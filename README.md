# dump2vtk（日本語）

LIGGGHTS/LAMMPS の dump を **VTK/VTU** に高速変換する単一ファイル **GUI + CLI** ツールです。  
粒子ダンプ、pair/local（force chain）ダンプ、ヘッダー置換を 1 本で扱えます。

- 主な特徴: 並列処理・チャンク処理、force chain の **ストリーミング読込**、自動連番収集（数字が末尾でなく**中間**でも検出）、**Louvain** によるコミュニティ検出と**コミュニティ集約**。
- 必要最小依存: `numpy`（GUI は `PySide6`、任意で `vtk` / `networkx` / `python-louvain`）。
- 英語版: [README_en.md](./README_en.md) ／ ライセンス: [LICENCE](./LICENCE)

---

## 1. Python の実行方法

### 1.1 動作環境
- Python **3.8 以上**（推奨: 3.10+）

### 1.2 インストール
```bash
# 最小構成
pip install numpy

# GUI を使う場合
pip install PySide6

# オプション機能（任意）
pip install vtk networkx python-louvain
```

### 1.3 起動（GUI）
```bash
python dump2vtk.py
```
- GUI が起動します。ドラッグ＆ドロップでファイル追加、詳細設定は **Advanced** から。

### 1.4 起動（CLI）
3 つのサブコマンドを使用できます。各 OS で共通に `python dump2vtk.py ...` としてください。

**(a) 粒子 Dump → VTK（VTK legacy / python-vtk）**
```bash
python dump2vtk.py lpp DUMP_FILES... \
  -o OUTROOT --format {ascii|binary} \
  --cpunum N --chunksize K --no-overwrite
```

**(b) ヘッダー置換（`ITEM: ENTRIES ...`）**
```bash
python dump2vtk.py rename \
  -H "ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz Fnx Fny Fnz Ftx Fty Ftz" \
  inputs*.dump --inplace
```

**(c) force chain → VTK/VTU（Louvain + 集約）**
```bash
python dump2vtk.py force forcechain-*.dump \
  --vtk-format {vtk|vtu} --encoding {ascii|binary} --keep-periodic \
  --resolution 1.0 --seed 42 --write-pointdata \
  --outdir out/ --nan-fill 0.0 \
  --cpunum N --chunksize K --no-overwrite
```

---

## 2. GUI の使い方

GUI は **3 タブ構成**です（下のスクリーンショット参照）。

### 2.1 Particles: Dump → VTK

![Particles タブ（連番自動収集と並列変換）](1.png)

1. 左のリストへ dump を **1 つドラッグ＆ドロップ**します。  
   同名プレフィクスの**連番ファイルを自動収集**します（数字の位置は末尾でも中間でも可）。
2. 右上の **出力先フォルダ**を必要に応じて指定します（空欄時は入力と同じ場所）。
3. **Advanced settings** では **ASCII/BINARY**，並列数，チャンクサイズ，`--no-overwrite`（既存スキップ）を設定します。  
   *補足*: writer backend は **legacy（純 Python）** と **vtk（python-vtk バインディング）** を切替可能です。
4. **Export VTK** を押します。進捗バーとログに  
   `Start: particle VTK (parallel)... Done.` と表示されれば成功です。

---

### 2.2 Header Replace: ヘッダー置換

![Header Replace タブ（ITEM: ENTRIES ... の一括置換）](2.png)

1. pair/local 形式の dump（例: `fc0.dump, fc2000.dump, ...`）を投入します。  
   WSL パス（`\\wsl.localhost\...`）も扱えます。
2. **Replacement header** に以下の 1 行を入力します。  
   ```text
   ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz Fnx Fny Fnz Ftx Fty Ftz
   ```
3. **Apply header replacement** を押します。各ファイルが `-renamed.dump` として複製されます。

---

### 2.3 Force chain: フォースネットワーク（VTK/VTU 出力）

![Force chain タブ（バッチ処理と Louvain 集約）](3.png)

1. `*-renamed.dump` を **1 つ選ぶ**と、同プレフィクスの**連番**を自動収集します。
2. **Export network** を押します。各ステップの **VTK/VTU** を出力します。  
   ログに `Start: force chain batch (parallel)... All done.` が表示されれば完了です。

---

## 3. VTK の出力例（ParaView）

出力ファイルの表示例です。

![VTKファイル（粒子）の表示例](particle.png)  
*図: 粒子（VTK POLYDATA）*

![VTKファイル（force chain）の表示例](force_chain.png)  
*図: フォースネットワーク（線分）*

![VTKファイル（Louvain 法による力のネットワーク）の表示例](louvain.png)  
*図: Louvain によるコミュニティ彩色*

---

## 4. I/O の要点

- **粒子ダンプ**: `x/y/z` と `xs/ys/zs` を自動判定し、必要に応じて **unscale**。`...x/..y/...z` や `f_*/c_*/v_*[1..3]` は **ベクトル**として扱い、その他は **スカラー**として書き出します。  
  legacy VTK では `POINTS`・`VERTICES`・`POINT_DATA` に加え、各ステップの **バウンディングボックス**（RECTILINEAR_GRID）も出力します。python-vtk バックエンドでも同様の幾何を生成します。
- **force chain**: `ENTRIES` に **12 必須列**（`x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz`）を想定。任意列はスカラー/ベクトルとして自動検出します。  
  エッジは `(id1, id2)` で生成し、`force = ||(fx,fy,fz)||`、`connectionLength` を計算します。**Louvain** で `community` / `intra_comm` を付与し、`--write-pointdata` でノード側 `node_community`・`node_degree`・`node_force_sum` を出力可能。さらに各属性の**コミュニティ mean/sum** を縁に展開します。

---

## 5. 連番ファイルの自動収集

- **prefix + digits + suffix + extension（.gz 対応）** の規則で、**1 ファイル投入だけで姉妹ファイル**を探索します。  
- 数字が **中間**にあるパターン（例: `forcechain-16000-renamed.dmp`）にも対応します。  
- ワイルドカード（例: `forcechain-*.dmp`）指定も可。重複は自動で除去します。

---

## 6. 謝辞 / Acknowledgments

- ダンプ処理および VTK 出力の設計は **Pizza.py / LPP（LIGGGHTS post-processing）** の考え方を参考にしています。  
- 例は **compaction_LIGGGHTS** を参照しました。作者とコミュニティに感謝します。

---

## 7. ライセンス / License

本リポジトリは **Pizza.py/LPP 由来のコード**を含むため、配布全体としては **GPL-2.0** に基づきます。  
同時に、本プロジェクトで新規に作成した部分は **MIT License** でも提供します（ただし GPL 由来部分を含む配布では **GPL-2.0** が優先）。

- SPDX: `GPL-2.0-only AND MIT`  
- 詳細は [`LICENCE`](./LICENCE) を参照してください。

---


`YOUR_ORG` をあなたの GitHub ネームスペースに置き換えてください。
