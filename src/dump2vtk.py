#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dump2vtk.py — Single-file GUI + CLI tool (PySide6) for LIGGGHTS dump workflows

・3つのツールを1つに統合（タブ分け）
  1) Dump → VTK (粒子)                 … 旧 lpp.py 相当（**並列化**・チャンク・スキップ生成対応）
  2) ヘッダー変換                      … 旧 dump_rename.py 相当（ENTRIES置換）
  3) Dump → VTK (force chain)         … 旧 dump2force.py 相当（VTU/VTK出力＋Louvain＋**コミュニティ集約**）

今回の主な改善点
----------------
- **force chain のストリーミング読込（増分読込）**を実装。巨大ファイルでも低メモリで処理。
- **force chain でも並列・チャンク処理**に対応（粒子タブと同じ UX）。
- **同名プレフィクスの連番自動選択を改良**。`forcechain-16000-renamed.dmp` のように「数字＋固定サフィックス」の
  パターンでも、1つ選択するだけで姉妹ファイルを自動収集（Dump→VTK(force chain) でも粒子と同様に働きます）。
- 純粋なバグ/ヒヤリ修正:
  * dy の差分計算で誤った列名を参照していたコード片を削除（`y2` を正しく使用）。
  * 連番展開 `expand_series` の堅牢化（数字が末尾でなく**中間**にあるケースをサポート）。
  * いくつかの I/O・型まわりを微修正（NaN/Inf 処理、型検査、Windows での mp 安定化など）。

・GUI: PySide6（**入力リスト領域のドラッグ＆ドロップ**、進捗バー、黒背景ログ、詳細設定の折り畳み）
・CUI: サブコマンド (lpp / rename / force) を用意。オプションでヘッドレス実行可能。
・ファイル選択: dump_xyzr-0, dump_xyzr-2000, dump_xyzr-36000 ... といった番号違いの姉妹ファイルを
  1つ選ぶだけで全部自動選択（参照ボタン/ドラッグ＆ドロップのどちらでも）。force chain でも同様。

依存: numpy （GUI使用時のみ PySide6）
オプション（任意）: vtk（--backend vtk 選択時/VTKライブラリ書き出し）
オプション（任意）: networkx, python-louvain（force chainでコミュニティ検出）

ライセンス: Pizza.py 由来部分は元のGPLを継承、それ以外はMITに準拠。
"""

from __future__ import annotations

import argparse
import concurrent.futures as _fut
import functools
import glob
import gzip
import io
import math
import multiprocessing as mp
import os
import re
import struct
import sys
from dataclasses import dataclass
from math import floor
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Callable

# ---------- Optional heavy deps ----------
try:
    import numpy as np
except Exception as e:
    print("This tool requires numpy. Please install it: pip install numpy", file=sys.stderr)
    raise

# VTK Python bindings are optional (only used when --backend vtk is selected)
_VTK_AVAILABLE = False
_vtk = None
def _try_import_vtk():
    global _VTK_AVAILABLE, _vtk
    if _VTK_AVAILABLE:
        return True
    try:
        import vtk as _vtk_mod  # type: ignore
        _vtk = _vtk_mod
        _VTK_AVAILABLE = True
    except Exception:
        _VTK_AVAILABLE = False
        _vtk = None
    return _VTK_AVAILABLE

# Optional for Louvain in force-network mode
try:
    import networkx as nx  # type: ignore
except Exception:
    nx = None  # type: ignore
try:
    from community import community_louvain  # type: ignore
except Exception:
    try:
        import community as community_louvain  # type: ignore
    except Exception:
        community_louvain = None  # type: ignore

# =============================================================================
#                              Shared utilities
# =============================================================================

def println(*args, file=sys.stdout, end="\n", quiet=False):
    """Unified print that can be redirected; respects 'quiet'."""
    if quiet:
        return
    print(*args, file=file, end=end)
    try:
        file.flush()
    except Exception:
        pass

def _split_series_stem(stem: str) -> Optional[Tuple[str, str, str]]:
    """
    stem を 「prefix + 数字 + suffix」に分解して返す。数字が**末尾でなく中間**にある場合も検出する。
    例:
      forcechain-16000-renamed -> ("forcechain-", "16000", "-renamed")
      dump_xyzr-2000           -> ("dump_xyzr-", "2000", "")
    見つからなければ None。
    """
    m = re.search(r"^(.*?)(\d+)(.*)$", stem)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3)

def expand_series(paths: List[str]) -> List[str]:
    """
    「dump_xyzr-0, dump_xyzr-2000, dump_xyzr-36000 ...」のような連番姉妹を自動展開。
    1つのファイル or パターンから、同フォルダ内で「同じ接頭語 + 数字 + 同じ後置 + 同じ拡張子(.gz含む)」を収集。
    ※ 数字が**末尾**だけでなく**中間**にあるケース（例: forcechain-16000-renamed.dmp）にも対応。
    """
    out: List[str] = []
    seen = set()
    for p in paths:
        p = os.path.abspath(p)
        # ワイルドカードは素直に展開
        if any(ch in p for ch in "*?[]"):
            for g in glob.glob(p):
                g = os.path.abspath(g)
                if os.path.isfile(g) and g not in seen:
                    out.append(g); seen.add(g)
            continue
        if os.path.isdir(p):
            continue
        d = os.path.dirname(p)
        b = os.path.basename(p)
        # .gz を一旦剥いで拡張子を見る
        gz = b.endswith(".gz")
        base_no_gz = b[:-3] if gz else b
        root, ext = os.path.splitext(base_no_gz)  # "stem" と ".dmp" など
        trip = _split_series_stem(root)
        if not trip:
            # 連番化できない場合はそのまま
            if p not in seen:
                out.append(p); seen.add(p)
            continue
        pre, num, post = trip
        cand: List[Tuple[int, str]] = []
        for name in os.listdir(d):
            full = os.path.join(d, name)
            if not os.path.isfile(full):
                continue
            gz2 = name.endswith(".gz")
            base2 = name[:-3] if gz2 else name
            stem2, ext2 = os.path.splitext(base2)
            if ext2 != ext:
                continue
            trip2 = _split_series_stem(stem2)
            if not trip2:
                continue
            pre2, num2, post2 = trip2
            if pre2 == pre and post2 == post:
                try:
                    cand.append((int(num2), os.path.abspath(full)))
                except ValueError:
                    continue
        cand.sort(key=lambda kv: kv[0])
        for _, c in cand:
            if c not in seen:
                out.append(c); seen.add(c)
    return out

def ensure_outdir_default_from_inputs(inputs: List[str]) -> str:
    """入力ファイルのフォルダを既定の出力先として返す（空ならカレント）"""
    for p in inputs:
        if os.path.isfile(p):
            return os.path.dirname(os.path.abspath(p))
    return os.getcwd()

# =============================================================================
#                          Minimal dump reader (particles)
# =============================================================================

class Snap:
    __slots__ = ("time","tselect","natoms","nselect","aselect","xlo","xhi","ylo","yhi","zlo","zhi","atoms")
    def __init__(self):
        self.time=0; self.tselect=0; self.natoms=0; self.nselect=0
        self.aselect=np.zeros(0,dtype=int)
        self.xlo=self.xhi=self.ylo=self.yhi=self.zlo=self.zhi=0.0
        self.atoms: Optional[np.ndarray]=None

class DumpReader:
    """
    LAMMPSスタイルdump（原子座標）を簡易に読み込む。
    xs/x/ys/y/zs/z を自動判定し、必要なら unscale する。
    names マップを構築し、atoms(N, ncol) を保持。
    """
    def __init__(self, filelist: List[str], log: Callable[[str],None]=lambda s: None, quiet=False):
        self.files = filelist
        self.snaps: List[Snap]=[]
        self.names: Dict[str,int]={}
        self.scale_original=-1
        self.quiet=quiet
        self.log = log

    def read_all(self):
        for fpath in self.files:
            opener = gzip.open if fpath.endswith(".gz") else open
            mode   = "rt"
            with opener(fpath, mode) as f:
                while True:
                    s = self._read_one(f)
                    if s is None:
                        break
                    self.snaps.append(s)
                    if not self.quiet:
                        self.log(f"{s.time} ")
        # sort & unique
        self.snaps.sort(key=lambda s: s.time)
        uniq: List[Snap]=[]; seen=set()
        for s in self.snaps:
            if s.time not in seen:
                uniq.append(s); seen.add(s.time)
        self.snaps = uniq
        # select all
        for s in self.snaps:
            s.tselect=1; s.nselect=s.natoms; s.aselect=np.ones(s.natoms, dtype=int)

        # unscale if needed
        if self.scale_original==1 and all(k in self.names for k in ("x","y","z")):
            for s in self.snaps:
                x,y,z = self.names["x"], self.names["y"], self.names["z"]
                xprd=s.xhi-s.xlo; yprd=s.yhi-s.ylo; zprd=s.zhi-s.zlo
                a=s.atoms
                a[:,x]=s.xlo + a[:,x]*xprd
                a[:,y]=s.ylo + a[:,y]*yprd
                a[:,z]=s.zlo + a[:,z]*zprd

    def _read_one(self, f) -> Optional[Snap]:
        try:
            if not f.readline():
                return None
            tline = f.readline()
            if not tline: return None
            time = int(tline.split()[0])
            if not f.readline(): return None  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().split()[0])

            s = Snap(); s.time=time; s.natoms=natoms; s.aselect=np.zeros(natoms, dtype=int)

            # BOX
            hdr = f.readline()
            xyz_line = f.readline().split(); s.xlo, s.xhi = float(xyz_line[0]), float(xyz_line[1])
            xyz_line = f.readline().split(); s.ylo, s.yhi = float(xyz_line[0]), float(xyz_line[1])
            xyz_line = f.readline().split(); s.zlo, s.zhi = float(xyz_line[0]), float(xyz_line[1])

            # ATOMS header
            item = f.readline()  # e.g., "ITEM: ATOMS id type x y z"
            if not self.names:
                self.scale_original=-1
                xflag=yflag=zflag=-1
                cols = item.split()[2:]
                for i,w in enumerate(cols):
                    if w in ("x","xu"): xflag=0; self.names["x"]=i
                    elif w in ("xs","xsu"): xflag=1; self.names["x"]=i
                    elif w in ("y","yu"): yflag=0; self.names["y"]=i
                    elif w in ("ys","ysu"): yflag=1; self.names["y"]=i
                    elif w in ("z","zu"): zflag=0; self.names["z"]=i
                    elif w in ("zs","zsu"): zflag=1; self.names["z"]=i
                    else: self.names[w]=i
                if (xflag,yflag,zflag)==(0,0,0): self.scale_original=0
                if (xflag,yflag,zflag)==(1,1,1): self.scale_original=1

            if natoms>0:
                first = f.readline().split(); ncol=len(first); rows=[first]
                for _ in range(natoms-1): rows.append(f.readline().split())
                data = np.array(rows, dtype=float)
            else:
                data = np.zeros((0, len(self.names)), dtype=float)
            s.atoms=data
            return s
        except Exception:
            return None

# =============================================================================
#                      Scalars/Vectors detection for particles
# =============================================================================

def find_scalars_vectors(names: Dict[str,int]) -> Tuple[Dict[str,int], Dict[str,int]]:
    """Original Pizza.py inspired heuristic to collect scalar and 3-comp vector columns."""
    vectors: Dict[str,int]={}
    scalars: Dict[str,int]={}
    # reverse mapping
    inv = {idx:name for name,idx in names.items()}
    if inv:
        for i in range(max(inv)+1):
            if i not in inv: inv[i]=""
    regvx=re.compile(r".*x$"); regvy=re.compile(r".*y$"); regvz=re.compile(r".*z$")
    regf=re.compile(r"f_.*\[\d+\]"); regc=re.compile(r"c_.*\[\d+\]"); regv=re.compile(r"v_.*\[\d+\]")
    i=0; maxi=max(inv) if inv else -1
    while i<=maxi:
        nm=inv[i]
        # ...x/...y/...z triplet
        if i+2<=maxi and regvx.match(inv[i]) and regvy.match(inv[i+1]) and regvz.match(inv[i+2]):
            pre = inv[i][:-1] if len(inv[i])>1 else inv[i]
            suf=(inv[i][-1], inv[i+1][-1], inv[i+2][-1])
            if set(suf)==set("xyz"):
                ix=i+suf.index("x"); iy=i+suf.index("y"); iz=i+suf.index("z")
                vectors[pre]=ix; i+=3; continue
        # f_/c_/v_[1..3]
        for rg in (regf,regc,regv):
            m=rg.match(inv[i])
            if m:
                name=inv[i]
                try: n=int(name.split("[")[1].split("]")[0])
                except: n=0
                base=name.split("[")[0]
                next1=f"{base}[{n+1}]"; next2=f"{base}[{n+2}]"
                if i+2<=maxi and inv[i+1]==next1 and inv[i+2]==next2:
                    vectors[base]=i; i+=3; break
                else:
                    scalars[base]=i; i+=1; break
        else:
            if nm!="": scalars[nm]=i
            i+=1
    if "x" not in vectors:
        print("Error: vector x y z must exist in dump.", file=sys.stderr)
        sys.exit(1)
    return scalars, vectors

# =============================================================================
#                         Particle → VTK writers
# =============================================================================

def _vtk_write_header_legacy(f, comment: str, dataset: str, binary=False):
    if binary:
        f.write(b"# vtk DataFile Version 2.0\n")
        f.write((comment+"\n").encode("ascii"))
        f.write(b"BINARY\n")
        f.write((f"DATASET {dataset}\n").encode("ascii"))
    else:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(comment+"\n")
        f.write("ASCII\n")
        f.write(f"DATASET {dataset}\n")

def _vtk_write_floats(f, arr: Iterable[float], binary=False):
    if binary:
        for v in arr: f.write(struct.pack(">f", float(v)))
        f.write(b"\n")
    else:
        buf=[]; c=0
        for v in arr:
            buf.append(f"{float(v):.9g}"); c+=1
            if c==9: f.write(" ".join(buf)+"\n"); buf=[]; c=0
        if buf: f.write(" ".join(buf)+"\n")

def _vtk_write_ints(f, arr: Iterable[int], binary=False):
    if binary:
        for v in arr: f.write(struct.pack(">i", int(v)))
        f.write(b"\n")
    else:
        buf=[]; c=0
        for v in arr:
            buf.append(str(int(v))); c+=1
            if c==16: f.write(" ".join(buf)+"\n"); buf=[]; c=0
        if buf: f.write(" ".join(buf)+"\n")

def write_particles_legacy_polydata(file: str,
                                    atoms_full: Optional[np.ndarray],
                                    names: Dict[str,int],
                                    fmt: str="ascii"):
    """Write particle points with scalars/vectors as POLYDATA (legacy VTK)."""
    binary = (fmt.lower()=="binary")
    npts = 0 if atoms_full is None else atoms_full.shape[0]
    scalars, vectors = find_scalars_vectors(names)

    mode = "wb" if binary else "w"
    with open(file, mode) as f:
        _vtk_write_header_legacy(f, "Generated by dump2vtk (legacy)", "POLYDATA", binary=binary)
        # POINTS
        if binary:
            f.write(f"POINTS {npts} float\n".encode("ascii"))
        else:
            f.write(f"POINTS {npts} float\n")
        coords=[]
        if npts:
            xidx = vectors["x"]
            coords = atoms_full[:, xidx:xidx+3].reshape(-1).tolist()
        _vtk_write_floats(f, coords, binary=binary)
        # VERTICES
        if binary:
            f.write(f"VERTICES {npts} {2*npts}\n".encode("ascii"))
        else:
            f.write(f"VERTICES {npts} {2*npts}\n")
        ints=[]
        for i in range(npts): ints.extend([1,i])
        _vtk_write_ints(f, ints, binary=binary)
        # POINT_DATA
        if binary: f.write(f"POINT_DATA {npts}\n".encode("ascii"))
        else:      f.write(f"POINT_DATA {npts}\n")
        # VECTORS (skip coords key 'x')
        for key, start in vectors.items():
            if key=="x": continue
            if binary: f.write(f"VECTORS {key} float\n".encode("ascii"))
            else:      f.write(f"VECTORS {key} float\n")
            vec = (atoms_full[:, start:start+3] if npts else np.zeros((0,3))).reshape(-1).tolist()
            _vtk_write_floats(f, vec, binary=binary)
        # SCALARS (int/float を自動判定)
        for key, idx in scalars.items():
            vals = (atoms_full[:, idx] if npts else np.array([], dtype=float))
            is_int = False
            if vals.size>0:
                # 「全要素が |v - round(v)| < 1e-12 」なら整数扱い
                dif = np.abs(vals - np.rint(vals))
                is_int = bool(np.all(dif < 1e-12))
            if binary:
                if is_int:
                    f.write(f"SCALARS {key} int 1\n".encode("ascii"))
                else:
                    f.write(f"SCALARS {key} float 1\n".encode("ascii"))
                f.write(b"LOOKUP_TABLE default\n")
            else:
                if is_int:
                    f.write(f"SCALARS {key} int 1\n")
                else:
                    f.write(f"SCALARS {key} float 1\n")
                f.write("LOOKUP_TABLE default\n")
            if is_int:
                _vtk_write_ints(f, [int(v) for v in vals.tolist()], binary=binary)
            else:
                _vtk_write_floats(f, vals.tolist(), binary=binary)

def write_bounding_box_legacy(file: str, xlo,xhi,ylo,yhi,zlo,zhi, fmt="ascii"):
    binary = (fmt.lower()=="binary")
    mode = "wb" if binary else "w"
    with open(file, mode) as f:
        _vtk_write_header_legacy(f, "Generated by dump2vtk (legacy)", "RECTILINEAR_GRID", binary=binary)
        line = "DIMENSIONS 2 2 2\n"
        f.write(line.encode("ascii") if binary else line)
        for axis,name in zip([(xlo,xhi),(ylo,yhi),(zlo,zhi)], ["X","Y","Z"]):
            line = f"{name}_COORDINATES 2 float\n"
            f.write(line.encode("ascii") if binary else line)
            _vtk_write_floats(f, axis, binary=binary)

def generate_filename(root: str, timestep: int) -> Tuple[str,str]:
    """Return (particles.vtk, bounding.vtk) using root prefix and timestep with zero-padding."""
    if   timestep<10:  s=f"{root}000{timestep}"
    elif timestep<100: s=f"{root}00{timestep}"
    elif timestep<1000:s=f"{root}0{timestep}"
    else:              s=f"{root}{timestep}"
    return s+".vtk", s+"_boundingBox.vtk"

# VTK backend using python-vtk (optional)
def write_particles_vtklib(file: str, atoms_full: Optional[np.ndarray], names: Dict[str,int], fmt="ascii"):
    if not _try_import_vtk():
        raise RuntimeError("vtk backend not available")
    vtk=_vtk
    npts = 0 if atoms_full is None else atoms_full.shape[0]
    poly = vtk.vtkPolyData(); points=vtk.vtkPoints(); verts=vtk.vtkCellArray()
    pdata = poly.GetPointData()
    if npts and "x" in names:
        xidx=names["x"]; coords = atoms_full[:, xidx:xidx+3]
        for i in range(npts):
            points.InsertNextPoint(float(coords[i,0]), float(coords[i,1]), float(coords[i,2]))
            verts.InsertNextCell(1); verts.InsertCellPoint(i)
    poly.SetPoints(points); poly.SetVerts(verts)
    scalars, vectors = find_scalars_vectors(names)
    # Vectors (skip 'x')
    for k, start in vectors.items():
        if k=="x": continue
        arr=vtk.vtkFloatArray(); arr.SetName(k); arr.SetNumberOfComponents(3)
        for i in range(npts):
            v=atoms_full[i, start:start+3]; arr.InsertNextTuple3(float(v[0]),float(v[1]),float(v[2]))
        pdata.AddArray(arr)
    # Scalars (int/float 切り替え)
    for k, idx in scalars.items():
        col = atoms_full[:, idx] if npts else np.array([], dtype=float)
        is_int = (col.size>0) and np.all(np.abs(col - np.rint(col)) < 1e-12)
        if is_int:
            arr=_vtk.vtkIntArray()
            arr.SetName(k)
            for v in col: arr.InsertNextValue(int(v))
        else:
            arr=vtk.vtkFloatArray(); arr.SetName(k)
            for v in col: arr.InsertNextValue(float(v))
        pdata.AddArray(arr)
    w=vtk.vtkPolyDataWriter(); w.SetFileName(file); w.SetInputData(poly)
    if fmt.lower()=="binary": w.SetFileTypeToBinary()
    else: w.SetFileTypeToASCII()
    w.Write()

def write_bb_vtklib(file: str, xlo,xhi,ylo,yhi,zlo,zhi, fmt="ascii"):
    if not _try_import_vtk():
        raise RuntimeError("vtk backend not available")
    vtk=_vtk
    grid=vtk.vtkRectilinearGrid(); grid.SetDimensions(2,2,2)
    def arr(vals):
        a=vtk.vtkFloatArray()
        for v in vals: a.InsertNextValue(float(v))
        return a
    grid.SetXCoordinates(arr([xlo,xhi])); grid.SetYCoordinates(arr([ylo,yhi])); grid.SetZCoordinates(arr([zlo,zhi]))
    w=vtk.vtkRectilinearGridWriter(); w.SetFileName(file); w.SetInputData(grid)
    if fmt.lower()=="binary": w.SetFileTypeToBinary()
    else: w.SetFileTypeToASCII()
    w.Write()

# ---------------- LPP (multiprocessing) ----------------

def _derive_output_root(infile: str, outroot: str) -> str:
    base = os.path.basename(infile)
    base = re.sub(r'\.gz$', '', base)
    base = re.sub(r'\.[^.]+$', '', base)
    m = re.match(r'^(.*?)([-_])\d+$', base)
    default_root = (m.group(1)+m.group(2)) if m else (base+"-")
    if outroot == "":
        return default_root
    elif outroot.endswith(("/", "\\")):
        return outroot + default_root
    else:
        return outroot


def _lpp_worker_onefile(args: Tuple[str, str, str, str, bool]) -> Tuple[str, List[int]]:
    fpath, outroot, fmt, backend, no_overwrite = args
    # Try to get timestep from file (2nd line)
    try:
        opener = gzip.open if fpath.endswith(".gz") else open
        with opener(fpath, "rt") as ff:
            ff.readline()
            tline = ff.readline()
            timestep_hint = int(tline.strip().split()[0])
    except Exception:
        timestep_hint = None

    # Output directory handling
    if outroot.endswith(("/", "\\")):
        outdir = outroot
    elif outroot:
        # If user gave a file-like root, treat as directory prefix
        outdir = outroot
        if not outdir.endswith(os.sep):
            outdir += os.sep
    else:
        outdir = ""  # let writer use working directory if not provided

    # read
    reader = DumpReader([fpath], quiet=True)
    reader.read_all()

    # Decide naming mode:
    #  A) 単一スナップ かつ 入力名が <prefix><sep><number> なら → {入力名}.vtk 方式
    #  B) それ以外 → 旧来の {prefix-}{zero-padded timestep}.vtk 方式
    base = os.path.basename(fpath)
    stem = re.sub(r'\.gz$', '', base)
    stem = re.sub(r'\.[^.]+$', '', stem)  # remove extension
    m = re.match(r'^(.*?)([-_]?)(\d+)$', stem)
    single_snap = (len(reader.snaps) == 1)
    use_input_stem = single_snap and (m is not None)

    if use_input_stem:
        # ex) dump_xyzr-2000 -> dump_xyzr-2000.vtk
        fp_out = (outdir or "") + stem + ".vtk"
        bb_out = (outdir or "") + stem + "_boundingBox.vtk"
        if no_overwrite and os.path.exists(fp_out) and os.path.exists(bb_out):
            return (fpath, [reader.snaps[0].time])
        # Write
        s = reader.snaps[0]
        if backend=="vtk":
            if not _try_import_vtk():
                backend="legacy"
            else:
                write_particles_vtklib(fp_out, s.atoms, reader.names, fmt=fmt)
                write_bb_vtklib(bb_out, s.xlo,s.xhi,s.ylo,s.yhi,s.zlo,s.zhi, fmt=fmt)
        if backend=="legacy":
            write_particles_legacy_polydata(fp_out, s.atoms, reader.names, fmt=fmt)
            write_bounding_box_legacy(bb_out, s.xlo,s.xhi,s.ylo,s.yhi,s.zlo,s.zhi, fmt=fmt)
        return (fpath, [s.time])

    # Fallback (multi-snapshot or non-numbered name): prefix + zero-padded timestep
    root = _derive_output_root(fpath, outroot)
    written = []
    for s in reader.snaps:
        fp_out, bb_out = generate_filename(root, s.time)
        if no_overwrite and os.path.exists(fp_out) and os.path.exists(bb_out):
            written.append(s.time); continue
        if backend=="vtk":
            if not _try_import_vtk():
                backend="legacy"
            else:
                write_particles_vtklib(fp_out, s.atoms, reader.names, fmt=fmt)
                write_bb_vtklib(bb_out, s.xlo,s.xhi,s.ylo,s.yhi,s.zlo,s.zhi, fmt=fmt)
        if backend=="legacy":
            write_particles_legacy_polydata(fp_out, s.atoms, reader.names, fmt=fmt)
            write_bounding_box_legacy(bb_out, s.xlo,s.xhi,s.ylo,s.yhi,s.zlo,s.zhi, fmt=fmt)
        written.append(s.time)
    return (fpath, written)


def run_lpp_parallel(filelist: List[str], output_root: str, fmt: str="ascii", backend: str="legacy",
                     cpunum: int=mp.cpu_count(), chunksize: int=8, no_overwrite: bool=False,
                     log: Callable[[str],None]=lambda s: None,
                     progress: Callable[[int,int],None]=lambda a,b: None):
    # expand patterns and series
    expanded: List[str]=[]
    for pat in filelist:
        if any(ch in pat for ch in "*?[]"):
            expanded.extend(glob.glob(pat))
        else:
            expanded.append(pat)
    expanded = expand_series(expanded)
    if not expanded:
        raise RuntimeError("no dump file specified")
    # chunking
    n = len(expanded)
    total = n
    done = 0
    log("Start: particle VTK (parallel)...\n")
    with mp.Pool(processes=max(1, cpunum)) as pool:
        args = [(f, output_root, fmt, backend, no_overwrite) for f in expanded]
        for _ in pool.imap_unordered(_lpp_worker_onefile, args, chunksize=max(1, chunksize)):
            done += 1
            progress(done, total)
    log("Done.\n")

# =============================================================================
#                         Header replacer (rename)
# =============================================================================

def replace_header(src_path: Path, header: str, inplace: bool, suffix: str, encoding: str) -> Path:
    pattern = re.compile(r"^ITEM: ENTRIES ")
    if inplace:
        with src_path.open("r", encoding=encoding, errors="replace") as fin, \
             io.open(str(src_path)+".tmp", "w", encoding=encoding) as tmp:
            for line in fin:
                if pattern.match(line): tmp.write(header+"\n")
                else: tmp.write(line)
        os.replace(str(src_path)+".tmp", str(src_path))
        return src_path
    else:
        dst = src_path.with_name(src_path.stem + suffix + src_path.suffix)
        with src_path.open("r", encoding=encoding, errors="replace") as fin, \
             dst.open("w", encoding=encoding) as fout:
            for line in fin:
                if pattern.match(line): fout.write(header+"\n")
                else: fout.write(line)
        return dst

# =============================================================================
#               Pair/Local dump → Force network (VTU/VTK) writers
# =============================================================================

REQUIRED12 = ["x1","y1","z1","x2","y2","z2","id1","id2","periodic","fx","fy","fz"]

@dataclass
class FSnap:
    time: int = 0
    natoms: int = 0
    atoms: Optional[np.ndarray] = None  # (natoms, ncol)
    colnames: List[str] = None

class BDump:
    """
    Pair/local dump reader supporting .gz **and** 増分読込（ストリーミング）。
    * read_all(): 互換の一括読み込み
    * next():     次のスナップショットだけ読み込む（.gz でも consumed 文字を管理）
    """
    def __init__(self, filename: str, read_all: bool=False):
        self.filename = filename
        self.snaps: List[FSnap]=[]
        self.increment = 0 if read_all else 1
        self.nextfile = 0
        self.eof = 0  # バイト or 文字単位の消費位置（gzは文字）
        if read_all:
            self.read_all()

    # ---- 増分 ----
    def next(self) -> int:
        if not self.increment:
            raise RuntimeError("BDump created in non-incremental mode")
        fname = self.filename
        is_gz = fname.endswith(".gz")
        if is_gz:
            with gzip.open(fname, "rt") as f:
                f.read(self.eof)  # 文字数ぶんスキップ
                snap, consumed = self._read_snapshot_text(f)
                if not snap:
                    return -1
                self.eof += consumed
        else:
            with open(fname, "rb") as fb:
                fb.seek(self.eof)
                snap = self._read_snapshot_binary_text(fb)
                if not snap:
                    return -1
                self.eof = fb.tell()
        # 重複タイムステップはスキップ
        if any(s.time == snap.time for s in self.snaps):
            return self.next()
        self.snaps.append(snap)
        return snap.time

    # ---- 一括 ----
    def read_all(self):
        opener = gzip.open if self.filename.endswith(".gz") else open
        mode = "rt" if self.filename.endswith(".gz") else "rb"
        with opener(self.filename, mode) as f:
            while True:
                s = self._read_snapshot_text(f) if mode=="rt" else self._read_snapshot_binary_text(f)
                if not s:
                    break
                if mode=="rt":
                    s = s[0]  # (snap, consumed)
                self.snaps.append(s)
        # ソート＋uniq
        self.snaps.sort(key=lambda s: s.time)
        uniq: List[FSnap]=[]
        for s in self.snaps:
            if not uniq or s.time != uniq[-1].time:
                uniq.append(s)
        self.snaps = uniq

    # ---- 内部：1スナップ読取り ----
    def _read_snapshot_text(self, f) -> Tuple[Optional[FSnap], int]:
        consumed=0
        def rline():
            nonlocal consumed
            s=f.readline(); consumed+=len(s); return s
        try:
            if not rline(): return None, consumed
            t=int(rline().split()[0])
            rline(); nat=int(rline().split()[0])
            header=rline(); 
            if header.startswith("ITEM: BOX BOUNDS"):
                rline(); rline(); rline()
                header=rline()
                if not header: return None, consumed
            toks=header.strip().split()
            try:
                idx=toks.index("ENTRIES"); col=toks[idx+1:]
            except: col=toks[2:]
            s=FSnap(); s.time=t; s.natoms=nat; s.colnames=col or []
            if nat>0:
                first=rline().split(); ncol=len(first)
                rows=[first]
                for _ in range(nat-1): rows.append(rline().split())
                data=np.array(rows, dtype=np.float64)
                s.atoms=data
            else:
                s.atoms=np.zeros((0, len(col)), dtype=np.float64)
            return s, consumed
        except Exception:
            return None, consumed

    def _read_snapshot_binary_text(self, fb) -> Optional[FSnap]:
        try:
            if not fb.readline(): return None
            t=int(fb.readline().split()[0])
            if not fb.readline(): return None
            nat=int(fb.readline().split()[0])
            header=fb.readline()
            if not header: return None
            try: header_s=header.decode("utf-8", errors="ignore").strip()
            except: header_s=str(header).strip()
            if header_s.startswith("ITEM: BOX BOUNDS"):
                fb.readline(); fb.readline(); fb.readline()
                header=fb.readline()
                if not header: return None
                try: header_s=header.decode("utf-8", errors="ignore").strip()
                except: header_s=str(header).strip()
            toks=header_s.split()
            try:
                idx=toks.index("ENTRIES"); col=toks[idx+1:]
            except: col=toks[2:]
            s=FSnap(); s.time=t; s.natoms=nat; s.colnames=col or []
            if nat>0:
                first=fb.readline().split(); ncol=len(first)
                rows=[first]; 
                for _ in range(nat-1): rows.append(fb.readline().split())
                data=np.array(rows, dtype=np.float64); s.atoms=data
            else:
                s.atoms=np.zeros((0,len(col)), dtype=np.float64)
            return s
        except Exception:
            return None

def _parse_required(colnames: List[str]) -> Dict[str,int]:
    m={}
    for r in REQUIRED12:
        if r not in colnames:
            raise RuntimeError(f"Required column '{r}' not found in ENTRIES header")
        m[r]=colnames.index(r)
    return m

def _split_optional(colnames: List[str]) -> Tuple[List[str], Dict[str,Tuple[int,int,int]]]:
    if len(colnames)<=12: return [], {}
    start=12; n=len(colnames); covered=set()
    scalars=[]; vectors={}
    def ends_xyz(name): 
        if not name: return None
        if name[-1].lower() in ("x","y","z"): return name[:-1]
        return None
    re_num=re.compile(r"^(?P<prefix>.+?)(?P<axis>[123])$")
    re_brk=re.compile(r"^(?P<prefix>.+?)\[(?P<idx>\d+)\]$")
    i=start
    while i<=n-3:
        a,b,c=colnames[i], colnames[i+1], colnames[i+2]
        pa,pb,pc = ends_xyz(a), ends_xyz(b), ends_xyz(c)
        if pa and pb and pc:
            suf=(a[-1].lower(), b[-1].lower(), c[-1].lower())
            if set(suf)=={"x","y","z"} and pa==pb==pc:
                ix=i+suf.index("x"); iy=i+suf.index("y"); iz=i+suf.index("z")
                vectors[pa]=(ix,iy,iz); covered.update({colnames[ix],colnames[iy],colnames[iz]}); i+=3; continue
        ma, nb, cc = re_num.match(a), re_num.match(b), re_num.match(c)
        if ma and nb and cc and ma.group("prefix")==nb.group("prefix")==cc.group("prefix") and {ma.group("axis"), nb.group("axis"), cc.group("axis")}=={"1","2","3"}:
            trip={ma.group("axis"):i, nb.group("axis"):i+1, cc.group("axis"):i+2}
            vectors[ma.group("prefix")]=(trip["1"],trip["2"],trip["3"]); covered.update({colnames[trip["1"]], colnames[trip["2"]], colnames[trip["3"]]}); i+=3; continue
        ma, nb, cc = re_brk.match(a), re_brk.match(b), re_brk.match(c)
        if ma and nb and cc and ma.group("prefix")==nb.group("prefix")==cc.group("prefix") and {ma.group("idx"), nb.group("idx"), cc.group("idx")}=={"1","2","3"}:
            trip={ma.group("idx"):i, nb.group("idx"):i+1, cc.group("idx"):i+2}
            ix,iy,iz = trip["1"], trip["2"], trip["3"]
            vectors[ma.group("prefix")]=(ix,iy,iz); covered.update({colnames[ix], colnames[iy], colnames[iz]}); i+=3; continue
        i+=1
    for k in range(start,n):
        name=colnames[k]
        if name not in covered: scalars.append(name)
    return scalars, vectors

def _build_force_geometry(s: FSnap, keep_periodic: bool):
    a = s.atoms if s.atoms is not None else np.zeros((0,0), dtype=np.float64)
    col = s.colnames or []
    idx = _parse_required(col)
    scalars_ex, vectors_ex = _split_optional(col)
    periodic = a[:, idx["periodic"]].astype(np.int64)!=0
    mask = np.ones(a.shape[0], dtype=bool) if keep_periodic else ~periodic
    id1_all = a[:, idx["id1"]].astype(np.int64); id2_all = a[:, idx["id2"]].astype(np.int64)
    id1 = id1_all[mask]; id2 = id2_all[mask]; M = id1.shape[0]
    ids = np.unique(np.concatenate([id1,id2])); N = ids.size
    id2idx = {int(pid):i for i,pid in enumerate(ids)}
    # points (choose pos from whichever side first seen)
    points = np.zeros((N,3), dtype=np.float64)
    first1={}; first2={}
    for i,pid in enumerate(id1_all):
        if int(pid) not in first1: first1[int(pid)]=i
    for i,pid in enumerate(id2_all):
        if int(pid) not in first2: first2[int(pid)]=i
    for pid, pidx in id2idx.items():
        if pid in first1:
            j=first1[pid]; points[pidx,0]=a[j,idx["x1"]]; points[pidx,1]=a[j,idx["y1"]]; points[pidx,2]=a[j,idx["z1"]]
        else:
            j=first2[pid]; points[pidx,0]=a[j,idx["x2"]]; points[pidx,1]=a[j,idx["y2"]]; points[pidx,2]=a[j,idx["z2"]]
    lines = np.empty((M,2), dtype=np.int32)
    for k in range(M):
        lines[k,0]=id2idx[int(id1[k])]; lines[k,1]=id2idx[int(id2[k])]
    # base cell data
    dx=a[:, idx["x1"]]-a[:, idx["x2"]]
    dy=a[:, idx["y1"]]-a[:, idx["y2"]]
    dz=a[:, idx["z1"]]-a[:, idx["z2"]]
    length=np.sqrt(dx*dx+dy*dy+dz*dz)[mask].astype(np.float64, copy=False)
    fx=a[:,idx["fx"]]; fy=a[:,idx["fy"]]; fz=a[:,idx["fz"]]; fmag=np.sqrt(fx*fx+fy*fy+fz*fz)[mask].astype(np.float64, copy=False)
    cell_scalars={"force": fmag, "connectionLength": length}
    for nm in scalars_ex:
        colv = a[:, col.index(nm)][mask].astype(np.float64, copy=False)
        cell_scalars[nm]=colv
    cell_vectors={}
    for pre,(ix,iy,iz) in vectors_ex.items():
        vx=a[:,ix][mask].astype(np.float64, copy=False); vy=a[:,iy][mask].astype(np.float64, copy=False); vz=a[:,iz][mask].astype(np.float64, copy=False)
        cell_vectors[pre]=np.stack([vx,vy,vz], axis=1)
    return points, lines, cell_scalars, cell_vectors, ids.astype(np.int64), (id1,id2)

def _louvain(ids, id_pairs, weight, resolution=1.0, seed=42, write_pointdata=False, comm_targets: Optional[Dict[str,np.ndarray]]=None):
    if community_louvain is None or nx is None:
        raise SystemExit("Louvain requires 'python-louvain' and 'networkx'. Install: pip install python-louvain networkx")
    id1,id2=id_pairs; M=weight.shape[0]
    G=nx.Graph(); G.add_nodes_from([int(x) for x in ids.tolist()])
    for a,b,w in zip(id1.tolist(), id2.tolist(), weight.tolist()):
        w=float(w); 
        if w>0 and a!=b:
            if G.has_edge(int(a),int(b)): G[int(a)][int(b)]["weight"]+=w
            else: G.add_edge(int(a),int(b),weight=w)
    part=community_louvain.best_partition(G, weight="weight", resolution=resolution, random_state=seed)
    node_comm = np.array([int(part.get(int(pid), -1)) for pid in ids], dtype=np.int64)
    comm_map={int(k):int(v) for k,v in part.items()}
    comm_e=np.empty(M, dtype=np.int64); intra=np.zeros(M, dtype=np.int32)
    for k in range(M):
        ca=comm_map.get(int(id1[k]),-1); cb=comm_map.get(int(id2[k]),-1)
        if ca==cb and ca!=-1: comm_e[k]=ca; intra[k]=1
        else: comm_e[k]=-1; intra[k]=0
    point_ann={}
    if write_pointdata:
        deg=dict(G.degree()); st=dict(G.degree(weight="weight"))
        node_degree=np.array([int(deg.get(int(pid),0)) for pid in ids], dtype=np.int64)
        node_force_sum=np.array([float(st.get(int(pid),0.0)) for pid in ids], dtype=np.float64)
        point_ann={"node_community": node_comm.astype(np.float64),
                   "node_degree": node_degree.astype(np.float64),
                   "node_force_sum": node_force_sum.astype(np.float64)}
    cell_ann={"community": comm_e.astype(np.float64), "intra_comm": intra.astype(np.float64)}
    # community aggregates
    if comm_targets:
        from collections import defaultdict
        edges_by_comm: Dict[int, List[int]] = defaultdict(list)
        for idx_e, cid in enumerate(comm_e):
            if cid != -1:
                edges_by_comm[int(cid)].append(idx_e)
        for name, arr in comm_targets.items():
            arr = np.asarray(arr, dtype=np.float64)
            out_mean = np.full(M, np.nan, dtype=np.float64)
            out_sum  = np.full(M, np.nan, dtype=np.float64)
            for cid, idxs in edges_by_comm.items():
                vals = arr[idxs]
                if len(vals):
                    out_sum[idxs]  = float(np.nansum(vals))
                    out_mean[idxs] = float(np.nanmean(vals))
            cell_ann[f"comm_mean_{name}"] = out_mean
            cell_ann[f"comm_sum_{name}"]  = out_sum
    return cell_ann, point_ann

def _finite_scalar(a: np.ndarray, fill: float) -> np.ndarray:
    b=np.asarray(a, dtype=np.float32).copy(order="C"); m=~np.isfinite(b); b[m]=fill; return b
def _finite_vec(v: np.ndarray, fill: float) -> np.ndarray:
    b=np.asarray(v, dtype=np.float32).reshape(-1,3).copy(order="C"); m=~np.isfinite(b); b[m]=fill; return b

def write_vtk_lines_ascii(path, points, lines, cell_scalars, cell_vectors, point_scalars=None, point_vectors=None, nan_fill=0.0):
    point_scalars=point_scalars or {}; point_vectors=point_vectors or {}
    npts=points.shape[0]; nlines=lines.shape[0]; total_ints=nlines*3
    with open(path,"w",encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n"); f.write("Force network\n"); f.write("ASCII\n"); f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {npts} float\n")
        for i in range(npts): f.write(f"{points[i,0]:.15g} {points[i,1]:.15g} {points[i,2]:.15g}\n")
        f.write(f"LINES {nlines} {total_ints}\n")
        for i in range(nlines): f.write(f"2 {int(lines[i,0])} {int(lines[i,1])}\n")
        if point_scalars or point_vectors:
            f.write(f"POINT_DATA {npts}\n")
            for name, arr in point_scalars.items():
                f.write(f"SCALARS {name} float 1\n"); f.write("LOOKUP_TABLE default\n")
                out=_finite_scalar(arr, nan_fill)
                for v in out: f.write(f"{float(v):.7g}\n")
            for name, vec in point_vectors.items():
                f.write(f"VECTORS {name} float\n")
                out=_finite_vec(vec, nan_fill)
                for r in out: f.write(f"{float(r[0]):.7g} {float(r[1]):.7g} {float(r[2]):.7g}\n")
        f.write(f"CELL_DATA {nlines}\n")
        for name, arr in cell_scalars.items():
            f.write(f"SCALARS {name} float 1\n"); f.write("LOOKUP_TABLE default\n")
            out=_finite_scalar(arr, nan_fill)
            for v in out: f.write(f"{float(v):.7g}\n")
        for name, vec in cell_vectors.items():
            f.write(f"VECTORS {name} float\n")
            out=_finite_vec(vec, nan_fill)
            for r in out: f.write(f"{float(r[0]):.7g} {float(r[1]):.7g} {float(r[2]):.7g}\n")

def write_vtk_lines_binary(path, points, lines, cell_scalars, cell_vectors, point_scalars=None, point_vectors=None, nan_fill=0.0):
    point_scalars=point_scalars or {}; point_vectors=point_vectors or {}
    npts=points.shape[0]; nlines=lines.shape[0]; total_ints=nlines*3
    with open(path,"wb") as f:
        f.write(b"# vtk DataFile Version 3.0\n"); f.write(b"Force network\n"); f.write(b"BINARY\n"); f.write(b"DATASET POLYDATA\n")
        f.write(f"POINTS {npts} float\n".encode("ascii"))
        if npts:
            fb = points.astype(np.float32).byteswap().tobytes(order="C")
            f.write(fb); f.write(b"\n")
        f.write(f"LINES {nlines} {total_ints}\n".encode("ascii"))
        if nlines:
            rec=np.empty((nlines,3), dtype=">i4"); rec[:,0]=2; rec[:,1:]=lines.astype(np.int32)
            f.write(rec.tobytes(order="C")); f.write(b"\n")
        if point_scalars or point_vectors:
            f.write(f"POINT_DATA {npts}\n".encode("ascii"))
            for name, arr in point_scalars.items():
                f.write(f"SCALARS {name} float 1\n".encode("ascii")); f.write(b"LOOKUP_TABLE default\n")
                out=_finite_scalar(arr, nan_fill); f.write(out.byteswap().tobytes(order="C")); f.write(b"\n")
            for name, vec in point_vectors.items():
                f.write(f"VECTORS {name} float\n".encode("ascii"))
                out=_finite_vec(vec, nan_fill); f.write(out.byteswap().tobytes(order="C")); f.write(b"\n")
        f.write(f"CELL_DATA {nlines}\n".encode("ascii"))
        for name, arr in cell_scalars.items():
            f.write(f"SCALARS {name} float 1\n".encode("ascii")); f.write(b"LOOKUP_TABLE default\n")
            out=_finite_scalar(arr, nan_fill); f.write(out.byteswap().tobytes(order="C")); f.write(b"\n")
        for name, vec in cell_vectors.items():
            f.write(f"VECTORS {name} float\n".encode("ascii"))
            out=_finite_vec(vec, nan_fill); f.write(out.byteswap().tobytes(order="C")); f.write(b"\n")

def _vtu_ascii_arrays(f, tag, scalars, vectors, nitems):
    f.write(f'      <{tag}>\n')
    for name, arr in scalars.items():
        f.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
        if nitems:
            flat=" ".join(str(float(x)) for x in np.asarray(arr, dtype=np.float64).ravel())
            f.write("          "+flat+"\n")
        f.write('        </DataArray>\n')
    for name, vec in vectors.items():
        f.write(f'        <DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">\n')
        if nitems:
            flat=" ".join(str(float(x)) for x in np.asarray(vec, dtype=np.float64).reshape(-1,3).ravel())
            f.write("          "+flat+"\n")
        f.write('        </DataArray>\n')
    f.write(f'      </{tag}>\n')

def write_vtu_ascii(path, points, lines, cell_scalars, cell_vectors, point_scalars=None, point_vectors=None):
    point_scalars=point_scalars or {}; point_vectors=point_vectors or {}
    npts=points.shape[0]; nc=lines.shape[0]
    connectivity=lines.astype(np.int32).ravel(order="C")
    offsets=((np.arange(nc, dtype=np.int32)+1)*2)
    types=np.full(nc,3,dtype=np.uint8)
    def arr_ascii(a):
        a=np.asarray(a)
        if a.dtype.kind in "iu": return " ".join(str(int(x)) for x in a.ravel())
        return " ".join(str(float(x)) for x in a.ravel())
    with open(path,"w",encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{npts}" NumberOfCells="{nc}">\n')
        _vtu_ascii_arrays(f, "PointData", point_scalars, point_vectors, npts)
        _vtu_ascii_arrays(f, "CellData", cell_scalars, cell_vectors, nc)
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        if npts: f.write("          "+arr_ascii(points)+"\n")
        f.write('        </DataArray>\n'); f.write('      </Points>\n')
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        if nc: f.write("          "+arr_ascii(connectivity)+"\n")
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        if nc: f.write("          "+arr_ascii(offsets)+"\n")
        f.write('        </DataArray>\n')
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        if nc: f.write("          "+arr_ascii(types)+"\n")
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')
        f.write('    </Piece>\n'); f.write('  </UnstructuredGrid>\n'); f.write('</VTKFile>\n')


def write_vtu_binary(path: str,
                     points: np.ndarray, lines: np.ndarray,
                     cell_scalars: Dict[str,np.ndarray], cell_vectors: Dict[str,np.ndarray],
                     point_scalars: Optional[Dict[str,np.ndarray]] = None,
                     point_vectors: Optional[Dict[str,np.ndarray]] = None):
    """Binary VTU (appended raw) — backslashを含むf-string式を回避"""
    point_scalars = point_scalars or {}
    point_vectors = point_vectors or {}
    npts, nc = points.shape[0], lines.shape[0]
    connectivity = lines.astype(np.int32).ravel(order="C")
    offsets = ((np.arange(nc, dtype=np.int32) + 1) * 2).astype(np.int32)
    types = np.full(nc, 3, dtype=np.uint8)
    blocks: List[Tuple[str, bytes, Dict[str,str], str]] = []
    blocks.append(("Points", points.astype("<f8").ravel(order="C").tobytes(), {"type":"Float64","NumberOfComponents":"3"}, "Points"))
    blocks.append(("Cells/connectivity", connectivity.astype("<i4").tobytes(), {"type":"Int32","Name":"connectivity"}, "Cells"))
    blocks.append(("Cells/offsets", offsets.astype("<i4").tobytes(), {"type":"Int32","Name":"offsets"}, "Cells"))
    blocks.append(("Cells/types", types.astype("|u1").tobytes(), {"type":"UInt8","Name":"types"}, "Cells"))
    for name, arr in cell_scalars.items():
        blocks.append((f"CellData/{name}", np.asarray(arr, dtype="<f8").tobytes(order="C"), {"type":"Float64","Name":name}, "CellData"))
    for name, vec in cell_vectors.items():
        blocks.append((f"CellData/{name}", np.asarray(vec, dtype="<f8").reshape(-1,3).tobytes(order="C"), {"type":"Float64","Name":name,"NumberOfComponents":"3"}, "CellData"))
    for name, arr in (point_scalars or {}).items():
        blocks.append((f"PointData/{name}", np.asarray(arr, dtype="<f8").tobytes(order="C"), {"type":"Float64","Name":name}, "PointData"))
    for name, vec in (point_vectors or {}).items():
        blocks.append((f"PointData/{name}", np.asarray(vec, dtype="<f8").reshape(-1,3).tobytes(order="C"), {"type":"Float64","Name":name,"NumberOfComponents":"3"}, "PointData"))
    offsets_in_appended: Dict[str,int] = {}
    sizes: Dict[str,int] = {}
    off = 0
    for key, data, _meta, _sec in blocks:
        offsets_in_appended[key] = off
        sizes[key] = len(data)
        off += 4 + len(data)
    def _attrs(meta: Dict[str, str]) -> str:
        return " ".join([f'{k}="{v}"' for k, v in meta.items()])
    with open(path, "wb") as f:
        f.write(b'<?xml version="1.0"?>\n')
        f.write(b'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n')
        f.write(b'  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{npts}" NumberOfCells="{nc}">\n'.encode("utf-8"))
        # PointData headers
        f.write(b'      <PointData>\n')
        for key, _data, meta, sec in blocks:
            if sec != "PointData":
                continue
            attrs = _attrs(meta)
            f.write(f'        <DataArray {attrs} format="appended" offset="{offsets_in_appended[key]}"/>\n'.encode("utf-8"))
        f.write(b'      </PointData>\n')
        # CellData headers
        f.write(b'      <CellData>\n')
        for key, _data, meta, sec in blocks:
            if sec != "CellData":
                continue
            attrs = _attrs(meta)
            f.write(f'        <DataArray {attrs} format="appended" offset="{offsets_in_appended[key]}"/>\n'.encode("utf-8"))
        f.write(b'      </CellData>\n')
        # Points
        f.write(b'      <Points>\n')
        f.write(f'        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{offsets_in_appended["Points"]}"/>\n'.encode("utf-8"))
        f.write(b'      </Points>\n')
        # Cells
        f.write(b'      <Cells>\n')
        for key in ("Cells/connectivity","Cells/offsets","Cells/types"):
            meta = {"Cells/connectivity":('Int32','connectivity'),
                    "Cells/offsets":('Int32','offsets'),
                    "Cells/types":('UInt8','types')}[key]
            f.write(f'        <DataArray type="{meta[0]}" Name="{meta[1]}" format="appended" offset="{offsets_in_appended[key]}"/>\n'.encode("utf-8"))
        f.write(b'      </Cells>\n')
        # AppendedData
        f.write(b'  <AppendedData encoding="raw">\n_')
        for key, data, _meta, _sec in blocks:
            f.write(struct.pack("<I", sizes[key]))
            f.write(data)
        f.write(b'\n  </AppendedData>\n')
        f.write(b'  </UnstructuredGrid>\n')
        f.write(b'</VTKFile>\n')

def run_force_convert(dumpfile: str, vtk_format: str="vtk", encoding: str="ascii",
                      keep_periodic: bool=False, resolution: float=1.0, seed: int=42,
                      write_pointdata: bool=False, outdir: Optional[str]=None,
                      nan_fill: float=0.0, skip_existing: bool=False,
                      log: Callable[[str],None]=lambda s: None,
                      progress: Callable[[int,int],None]=lambda a,b: None):
    """1ファイルを**ストリーミング**でVTK/VTUに変換（Louvain付与）"""
    if not os.path.isfile(dumpfile):
        raise FileNotFoundError(dumpfile)
    base = os.path.basename(dumpfile)
    # 拡張子を剥がして「素の stem」を取得（.gz 二重拡張にも対応）
    stem = re.sub(r'\.gz$', '', base)
    stem = re.sub(r'\.[^.]+$', '', stem)  # 最後の拡張子 (.dmp 等) を落とす
    # もし 'dump.' などの接頭語を落としたい場合は次行を有効化
    # stem = re.sub(r'(?i)^dump[._-]', '', stem)
    outdir = outdir or os.path.dirname(os.path.abspath(dumpfile))
    if not os.path.isdir(outdir):
        raise RuntimeError(f"Output directory does not exist: {outdir}")
    b=BDump(dumpfile, read_all=False)
    total=0; written=0; fileindex=0
    log(f"Start: force chain convert ({os.path.basename(dumpfile)}) ...\n")
    # 1ステップずつ読んで処理
    while True:
        try:
            timestep = b.next()
        except Exception:
            break
        if timestep < 0:
            break
        total += 1
        snap = b.snaps[fileindex]
        # ファイル名中に既に数字がある場合は『その数字を現在の timestep に置換』
        trip = _split_series_stem(stem)  # -> (pre, num, post) / 数字が無ければ None
        if trip:
            pre, num, post = trip
            # 元の桁数を維持（例: 00010 → 00010）
            tstr = f"{int(timestep):0{len(num)}d}"
            vt_base = os.path.join(outdir, f"{pre}{tstr}{post}")
        else:
            # 数字が無い stem のときだけ _{timestep} を付ける
            vt_base = os.path.join(outdir, f"{stem}_{timestep}")
        vt_path = vt_base + (".vtk" if vtk_format=="vtk" else ".vtu")
        if skip_existing and os.path.exists(vt_path):
            progress(total, total)  # 保守: 1ファイルにつき1ステップのみのケースでも100%に
            fileindex += 1
            continue
        # 幾何＋基本属性
        points, lines, cell_scalars, cell_vectors, ids, id_pairs = _build_force_geometry(snap, keep_periodic=keep_periodic)
        # Louvain + community aggregates
        try:
            comm_targets = {**cell_scalars, **{k: np.linalg.norm(v,axis=1) for k,v in cell_vectors.items()}}
            louv_cell, louv_point = _louvain(ids, id_pairs, cell_scalars["force"], resolution=resolution, seed=seed, write_pointdata=write_pointdata, comm_targets=comm_targets)
        except SystemExit as e:
            raise RuntimeError(str(e))
        for k,v in louv_cell.items(): cell_scalars[k]=v
        point_scalars = louv_point if write_pointdata else {}
        # 書き出し
        if vtk_format=="vtu":
            if encoding=="ascii":
                write_vtu_ascii(vt_base+".vtu", points, lines, cell_scalars, cell_vectors, point_scalars, {})
            else:
                write_vtu_binary(vt_base+".vtu", points, lines, cell_scalars, cell_vectors, point_scalars, {})
        else:
            if encoding=="ascii":
                write_vtk_lines_ascii(vt_base+".vtk", points, lines, cell_scalars, cell_vectors, point_scalars, {}, nan_fill=nan_fill)
            else:
                write_vtk_lines_binary(vt_base+".vtk", points, lines, cell_scalars, cell_vectors, point_scalars, {}, nan_fill=nan_fill)
        written+=1
        log(f"t={timestep}: points={points.shape[0]} lines={lines.shape[0]} -> {os.path.basename(vt_path)}\n")
        progress(written, total)  # ここではファイル内スナップ数が未知のため、相対進捗として扱う
        fileindex += 1
    log(f"{dumpfile}: wrote {written}/{total}\n")

def _force_worker(args):
    (infile, params) = args
    try:
        run_force_convert(infile, **params, log=lambda s: None, progress=lambda a,b: None)
        return (infile, True, "")
    except Exception as e:
        return (infile, False, str(e))

def run_force_batch(files: List[str],
                    vtk_format: str="vtk", encoding: str="ascii",
                    keep_periodic: bool=False, resolution: float=1.0, seed: int=42,
                    write_pointdata: bool=False, outdir: Optional[str]=None,
                    nan_fill: float=0.0, skip_existing: bool=False, cpunum: int=mp.cpu_count(),
                    chunksize: int=8,
                    log: Callable[[str],None]=lambda s: None,
                    progress: Callable[[int,int],None]=lambda a,b: None):
    """force chain の複数ファイルを**並列**で処理。imap_unordered の chunksize を指定可能。"""
    # expand and uniq
    expanded: List[str]=[]
    for pat in files:
        if any(ch in pat for ch in "*?[]"):
            expanded.extend(glob.glob(pat))
        else:
            expanded.append(pat)
    expanded = expand_series(expanded)  # ★ 連番も自動収集（数字が中間でもOK）
    expanded = sorted(set(expanded))
    if not expanded:
        raise RuntimeError("no input file(s)")
    total=len(expanded); done=0
    log("Start: force chain batch (parallel)...\n")
    params = dict(vtk_format=vtk_format, encoding=encoding, keep_periodic=keep_periodic, resolution=resolution, seed=seed,
                  write_pointdata=write_pointdata, outdir=outdir, nan_fill=nan_fill, skip_existing=skip_existing)
    with mp.Pool(processes=max(1, cpunum)) as pool:
        for infile, ok, err in pool.imap_unordered(_force_worker, [(f, params) for f in expanded], chunksize=max(1, int(chunksize))):
            done += 1
            if ok: log(f"✔ {os.path.basename(infile)} done\n")
            else:  log(f"✖ {os.path.basename(infile)} failed: {err}\n")
            progress(done, total)
    log("All done.\n")

# =============================================================================
#                                   CLI
# =============================================================================

def build_cli():
    p=argparse.ArgumentParser(description="dump2vtk: GUI+CLI unified tool for LIGGGHTS dump conversion")
    sub=p.add_subparsers(dest="cmd", help="If omitted, GUI will be launched.")
    # lpp
    lp=sub.add_parser("lpp", help="Dump (particles) → VTK/BoundingBox output (VTK legacy/VTKlib)")
    lp.add_argument("dump_files", nargs="+", help="dump files or pattern (series auto-expansion)")
    lp.add_argument("-o","--output", dest="output_root", default="", help="Output file name/prefix (default: derived from first input)")
    lp.add_argument("--format", choices=["ascii","binary"], default="ascii", help="VTK legacy encoding")
    lp.add_argument("--backend", choices=["legacy","vtk"], default="legacy", help="writer backend ('vtk' requires python-vtk)")
    lp.add_argument("--cpunum", type=int, default=max(1, mp.cpu_count()), help="number of processes")
    lp.add_argument("--chunksize", type=int, default=8, help="imap_unordered chunk size")
    lp.add_argument("--no-overwrite", action="store_true", help="skip if .vtk already exists")
    # rename
    rn=sub.add_parser("rename", help="Replace header 'ITEM: ENTRIES ...'")
    rn.add_argument("-H","--header", required=True, help="Replacement header (e.g. 'ITEM: ENTRIES x1 y1 ...')")
    rn.add_argument("inputs", nargs="+", help="Input files/patterns")
    rn.add_argument("--inplace", action="store_true", help="Overwrite in-place. If not set, write to a new file with suffix '-renamed'")
    rn.add_argument("--suffix", default="-renamed", help="Suffix when not in-place (default: -renamed)")
    rn.add_argument("--encoding", default="utf-8", help="Text encoding (default: utf-8)")
    # force
    fc=sub.add_parser("force", help="Dump (pair/local) → line network (VTU/VTK) with Louvain (parallel multi-file)")
    fc.add_argument("dump_files", nargs="+", help="Input dump(s) (pick one to auto-collect numbered siblings)")
    fc.add_argument("--vtk-format", choices=["vtu","vtk"], default="vtk", help="Output format (default: vtk)")
    fc.add_argument("--encoding", choices=["ascii","binary"], default="ascii", help="VTK ASCII/BINARY")
    fc.add_argument("--keep-periodic", action="store_true", help="Keep periodic contacts")
    fc.add_argument("--resolution", type=float, default=1.0, help="Louvain resolution")
    fc.add_argument("--seed", type=int, default=42, help="Louvain random seed")
    fc.add_argument("--write-pointdata", action="store_true", help="Also write PointData (node_community/degree/force_sum)")
    fc.add_argument("--outdir", default=None, help="Output directory (default: same as input)")
    fc.add_argument("--nan-fill", type=float, default=0.0, help="(VTK legacy) Replace NaN/Inf with this value (default 0.0)")
    fc.add_argument("--cpunum", type=int, default=max(1, mp.cpu_count()), help="number of processes")
    fc.add_argument("--chunksize", type=int, default=8, help="imap_unordered chunk size (default 8)")
    fc.add_argument("--no-overwrite", action="store_true", help="Skip existing files (per timestep)")
    fc.add_argument("--quiet", action="store_true", help="Suppress progress/log outputs")
    return p

def run_cli(args: argparse.Namespace):
    if args.cmd=="lpp":
        def log(s): print(s, end="")
        def progress(i, n):
            pct = 0 if n==0 else int(i*100/n)
            print(f"\r[{pct:3d}%] ", end="")
            if i==n: print("")
        run_lpp_parallel(args.dump_files, args.output_root, fmt=args.format, backend=args.backend,
                         cpunum=args.cpunum, chunksize=args.chunksize, no_overwrite=args.no_overwrite,
                         log=log, progress=progress)
        return 0
    elif args.cmd=="rename":
        files=[]
        for pat in args.inputs:
            got = glob.glob(pat)
            if got: files.extend(got)
            elif os.path.exists(pat): files.append(pat)
            else: print(f"⚠️  No match: {pat}", file=sys.stderr)
        if not files:
            print("No input files found.", file=sys.stderr); return 1
        files=sorted(set(files))
        for f in files:
            try:
                out = replace_header(Path(f), args.header, args.inplace, args.suffix, args.encoding)
                print(f"✅ {f} -> {out}")
            except Exception as e:
                print(f"❌ Error processing {f}: {e}", file=sys.stderr)
        return 0
    elif args.cmd=="force":
        def log(s):
            if not args.quiet:
                print(s, end="")
        def progress(i,n):
            if not args.quiet:
                pct = 0 if n==0 else int(i*100/n)
                print(f"\r[{pct:3d}%] ", end="")
                if i==n: print("")
        run_force_batch(args.dump_files, vtk_format=args.vtk_format, encoding=args.encoding,
                        keep_periodic=args.keep_periodic, resolution=args.resolution, seed=args.seed,
                        write_pointdata=args.write_pointdata, outdir=args.outdir, nan_fill=args.nan_fill,
                        skip_existing=args.no_overwrite, cpunum=args.cpunum, chunksize=args.chunksize,
                        log=log, progress=progress)
        return 0
    else:
        # GUI fallback
        return run_gui()

# =============================================================================
#                                   GUI
# =============================================================================

def run_gui():
    """
    GUI (PySide6). Fixes & customizations:
      - 詳細設定の**日本語化**と**ツールチップ**。
      - 「並列数/チャンクサイズ」を**環境に応じて自動推奨**（ファイル選択時に再計算）。
      - 指定アイコン（input.png / clear.png / output.png / setting.png / convert.png）に対応。
      - プルダウン選択時に文字が見えなくなる問題を**スタイルで修正**。
      - AA_UseHighDpiPixmaps の代替として setHighDpiScaleFactorRoundingPolicy をアプリ生成前に設定。
      - ドロップ領域のフォントサイズ安全化。
      - ★ 重要: 「詳細設定」トグルに setting.png を必ず使用し、かつ**カーソルを近づけると説明が即表示**されるよう eventFilter と QToolTip を追加。
    """
    try:
        from PySide6 import QtCore, QtGui, QtWidgets
    except Exception as e:
        print("PySide6 is required for GUI. Install: pip install PySide6", file=sys.stderr)
        return 1

    # ---------------- i18n (GUI only) ----------------
    class I18N(QtCore.QObject):
        langChanged = QtCore.Signal(str)
        def __init__(self, lang="en"):
            super().__init__()
            self.lang = "en" if lang not in ("ja","en") else lang
        def set_lang(self, lang: str):
            lang = "ja" if str(lang).lower().startswith("ja") else "en"
            if lang != self.lang:
                self.lang = lang
                self.langChanged.emit(self.lang)
        def tr(self, en: str, ja: str) -> str:
            return en if self.lang == "en" else ja

    i18n = I18N("en")  # default English

    # ---------------- HiDPI policy (call BEFORE app instance) ----------------
    try:
        QtGui.QGuiApplication.setHighDpiScaleFactorRoundingPolicy(
            QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
        )
    except Exception:
        pass

    # ---------------- helper: auto-parallel guess ----------------
    def guess_workers() -> int:
        n = mp.cpu_count() or 1
        if n >= 16: return n - 2
        if n >= 8:  return n - 1
        return max(1, n)
    def guess_chunksize(nfiles: int, workers: int) -> int:
        if workers <= 0: workers = 1
        if nfiles <= 0:
            return 8 if workers <= 8 else 12
        # ざっくり: (総ジョブ/ワーカー)/2 を上限 32 で丸め
        g = max(1, int(math.ceil(nfiles / (workers * 2.0))))
        return int(min(32, g))

    # ---------------- icons ----------------
    def icons_base_dir():
        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            return sys._MEIPASS
        return os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = icons_base_dir()
    # 追加: どこに展開されても見つけられるようにする
    def resource_path(*parts):
        base = os.path.dirname(os.path.abspath(__file__))  # Nuitka onefile でもここが展開先になる
        cand = os.path.join(base, *parts)
        if os.path.exists(cand):
            return cand
        # assets/ に入れている場合のフォールバック
        cand_assets = os.path.join(base, "assets", *parts)
        return cand_assets

    # 既存の load_icon を assets 対応に
    def load_icon(fname: str, fallback: Optional[QtGui.QIcon]=None) -> QtGui.QIcon:
        for p in (resource_path(fname),):
            if os.path.exists(p):
                return QtGui.QIcon(p)
        return fallback if fallback is not None else QtGui.QIcon()

    # ---------------- Style tokens (compact & modern) ----------------
    ACCENT="#2563eb"  # blue
    OK="#16a34a"      # green
    WARN="#f59e0b"    # amber
    ERR="#ef4444"     # red
    BG_DARK="#111315" # log bg
    FG_LIGHT="#e5e7eb"

    GLOBAL_QSS = f"""
    QWidget {{
        font-family: "Segoe UI", "Hiragino Sans", "Noto Sans CJK JP", "Meiryo", sans-serif;
        font-size: 12px;
        color: #111111;
    }}
    QTabWidget::pane {{
        border: 1px solid #e5e7eb; border-radius: 6px; background: #ffffff;
    }}
    QTabBar::tab {{
        padding: 6px 10px; border: 1px solid transparent; margin: 2px; border-radius: 4px;
    }}
    QTabBar::tab:selected {{
        background: #eef2ff; border-color: #c7d2fe;
    }}
    QLabel[hint='true'] {{
        color: #6b7280;
    }}
    QToolButton, QPushButton {{
        padding: 6px 10px; border-radius: 6px; border: 1px solid #e5e7eb; background: #ffffff;
    }}
    QPushButton[accent='true'] {{
        background: {ACCENT}; color: white; border: none; padding: 8px 14px;
    }}
    QPushButton[accent='true']:disabled {{
        background: #9ca3af;
    }}
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
        padding: 5px 7px; border: 1px solid #e5e7eb; border-radius: 6px; background: #ffffff; color: #111111;
    }}
    /* ▼ ドロップダウン（プルダウン）で選択文字が見えなくなる問題を回避 */
    QComboBox QAbstractItemView {{
        background: #ffffff; 
        border: 1px solid #e5e7eb; 
        outline: 0;
    }}
    QComboBox QAbstractItemView::item {{
        padding: 4px 8px; color: #111111;
    }}
    QComboBox QAbstractItemView::item:selected {{
        background: {ACCENT}; color: #ffffff;
    }}
    QListWidget {{
        border: 1px dashed #c7d2fe; border-radius: 8px; background: #fafafa;
    }}
    QProgressBar {{
        border: 1px solid #e5e7eb; border-radius: 6px; background: #f3f4f6; text-align: center;
        min-height: 16px;
    }}
    QProgressBar::chunk {{
        background-color: {OK}; border-radius: 6px;
    }}
    QGroupBox {{
        border: 1px solid #e5e7eb; border-radius: 8px; margin-top: 10px; padding: 10px;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin; left: 8px; padding: 0px 4px;
    }}
    /* ★ ツールチップの視認性を強化（見えない/出ない問題への対処） */
    QToolTip {{
        background-color: #111827;
        color: #f9fafb;
        border: 1px solid #374151;
        padding: 6px 8px;
        border-radius: 6px;
    }}
    """

    # ------------- Reusable widgets -------------
    class LogView(QtWidgets.QTextEdit):
        def __init__(self, parent=None):
            super().__init__(parent)
            self.setReadOnly(True)
            self.setAcceptDrops(False)
            self.setMinimumHeight(140)
            self.setStyleSheet(f"""
                QTextEdit {{
                    background: {BG_DARK};
                    color: {FG_LIGHT};
                    font-family: Consolas, 'SF Mono', Menlo, 'Courier New', monospace;
                    border: 1px solid #222;
                    padding: 8px;
                    border-radius: 8px;
                }}
            """)
        def append_colored(self, text: str, color: str="#cbd5e1"):
            cursor=self.textCursor()
            fmt=QtGui.QTextCharFormat(); fmt.setForeground(QtGui.QBrush(QtGui.QColor(color)))
            cursor.movePosition(QtGui.QTextCursor.End)
            cursor.insertText(text, fmt)
            self.setTextCursor(cursor); self.ensureCursorVisible()

    class DropList(QtWidgets.QListWidget):
        filesDropped = QtCore.Signal(list)
        def __init__(self, parent=None):
            super().__init__(parent)
            self.setAcceptDrops(True)
            self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
            self.setAlternatingRowColors(True)
            self._dragActive=False
            self._placeholder="Drag & Drop files here (or use the folder icon)"
        def setPlaceholder(self, text: str):
            self._placeholder = text
            self.viewport().update()
        def dragEnterEvent(self, e):
            if e.mimeData().hasUrls():
                e.acceptProposedAction()
                self._dragActive=True; self.viewport().update()
            else:
                e.ignore()
        def dragMoveEvent(self, e):
            if e.mimeData().hasUrls():
                e.acceptProposedAction()
            else:
                e.ignore()
        def dragLeaveEvent(self, e):
            self._dragActive=False; self.viewport().update(); e.accept()
        def dropEvent(self, e):
            self._dragActive=False; self.viewport().update()
            paths=[]
            for u in e.mimeData().urls():
                p=u.toLocalFile()
                if os.path.isfile(p): paths.append(p)
            if paths:
                self.filesDropped.emit(paths)
            e.acceptProposedAction()
        def paintEvent(self, ev):
            super().paintEvent(ev)
            if self.count()==0:
                p=QtGui.QPainter(self.viewport())
                p.setRenderHint(QtGui.QPainter.Antialiasing)
                rect=self.viewport().rect().adjusted(6,6,-6,-6)
                pen=QtGui.QPen(QtGui.QColor("#93c5fd" if self._dragActive else "#cbd5e1"))
                pen.setStyle(QtCore.Qt.DashLine)
                p.setPen(pen)
                p.drawRoundedRect(rect, 8,8)
                # ---- Guard against negative font point size ----
                font=p.font()
                ps = font.pointSizeF()
                if ps is not None and ps > 0:
                    new_ps = max(1.0, ps*0.95)
                    font.setPointSizeF(new_ps)
                else:
                    px = font.pixelSize()
                    if px is not None and px > 0:
                        font.setPixelSize(max(1, int(px*0.95)))
                    else:
                        font.setPointSize(11)  # safe default
                p.setFont(font)
                p.setPen(QtGui.QColor("#6b7280"))
                p.drawText(rect, QtCore.Qt.AlignCenter, self._placeholder)

    class Collapsible(QtWidgets.QWidget):
        def __init__(self, title="Advanced", parent=None):
            super().__init__(parent)
            self.toggle=QtWidgets.QToolButton(text=title, checkable=True, checked=False)
            self.toggle.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
            self.toggle.setArrowType(QtCore.Qt.RightArrow)
            self.toggle.setStyleSheet("QToolButton { font-weight: 600; }")
            # アイコン（setting.png）— 明示指定
            self.toggle.setIcon(load_icon("setting.png"))
            # ★ ツールチップが「発現しない」環境対策: 即時表示 + 長めの表示時間 + Hover有効化
            try:
                self.toggle.setToolTipDuration(10000)  # 10s
            except Exception:
                pass
            self.toggle.setMouseTracking(True)
            try:
                self.toggle.setAttribute(QtCore.Qt.WA_Hover, True)
            except Exception:
                pass
            self.toggle.installEventFilter(self)

            self.toggle.toggled.connect(self._on_tog)
            self.body=QtWidgets.QWidget()
            self.body.setVisible(False)
            lay=QtWidgets.QVBoxLayout(self); lay.setContentsMargins(0,0,0,0); lay.setSpacing(6)
            lay.addWidget(self.toggle); lay.addWidget(self.body)

        def eventFilter(self, obj, ev):
            # ホバーした瞬間に QToolTip を明示表示（OS/テーマに依存せず確実に表示）
            if obj is self.toggle and ev.type() == QtCore.QEvent.Enter:
                QtWidgets.QToolTip.showText(QtGui.QCursor.pos(), self.toggle.toolTip(), self.toggle)
            return super().eventFilter(obj, ev)

        def _on_tog(self, on):
            self.toggle.setArrowType(QtCore.Qt.DownArrow if on else QtCore.Qt.RightArrow)
            self.body.setVisible(on)

        def setContentLayout(self, l):
            self.body.setLayout(l)

        def setTitle(self, title: str):
            self.toggle.setText(title)

        def setToolTip(self, tip: str):
            # Collapsible 全体ではなくトグルボタンにツールチップを設定
            self.toggle.setToolTip(tip)

    class Runner(QtCore.QThread):
        progressed = QtCore.Signal(int,int)
        logged = QtCore.Signal(str, str)  # text, color
        finished_ok = QtCore.Signal()
        failed = QtCore.Signal(str)
        def __init__(self, fn: Callable, args: tuple, kwargs: dict):
            super().__init__()
            self.fn=fn; self.args=args; self.kwargs=kwargs
        def run(self):
            def log(txt, color="#cbd5e1"): self.logged.emit(txt, color)
            def prog(i,n): self.progressed.emit(i,n)
            try:
                self.kwargs["log"]=log
                self.kwargs["progress"]=prog
                self.fn(*self.args, **self.kwargs)
                self.finished_ok.emit()
            except Exception as e:
                self.failed.emit(str(e))

    class BaseTab(QtWidgets.QWidget):
        requestOutdirSet = QtCore.Signal(str)
        def __init__(self, bg_color: str, parent=None):
            super().__init__(parent)
            self.setAutoFillBackground(True)
            pal=self.palette(); pal.setColor(QtGui.QPalette.Window, QtGui.QColor(bg_color)); self.setPalette(pal)
            self.log = LogView()
            self.progress = QtWidgets.QProgressBar()
            self.progress.setRange(0,100); self.progress.setValue(0); self.progress.setTextVisible(True)
            self.progress.setFormat("%p% in progress")
            self.files = DropList()
            self.files.filesDropped.connect(self._on_files_dropped)
        def _on_files_dropped(self, paths):
            expanded = expand_series(paths)
            self.add_files(expanded)
            out = ensure_outdir_default_from_inputs(expanded)
            self.requestOutdirSet.emit(out)
        def add_files(self, paths: List[str]):
            for p in paths:
                if not any(self.files.item(i).text()==p for i in range(self.files.count())):
                    self.files.addItem(p)
        def take_all_files(self) -> List[str]:
            return [self.files.item(i).text() for i in range(self.files.count())]
        def clear_files(self):
            self.files.clear()
        def gui_log(self, text:str, color:str="#cbd5e1"):
            self.log.append_colored(text, color)
        def set_progress(self, i:int, n:int):
            val = 0 if n==0 else int(i*100/n)
            self.progress.setValue(val)

        # Compact header bar builder
        def build_header_bar(self, title: str, widgets_right: List[QtWidgets.QWidget]=None):
            bar = QtWidgets.QWidget()
            lay = QtWidgets.QHBoxLayout(bar); lay.setContentsMargins(0,0,0,0); lay.setSpacing(6)
            lab = QtWidgets.QLabel(title); lab.setProperty("hint","true")
            lay.addWidget(lab); lay.addStretch(1)
            if widgets_right:
                for w in widgets_right:
                    lay.addWidget(w)
            # return both bar and label for live i18n
            bar._title_label = lab  # type: ignore[attr-defined]
            return bar, lab

        def two_pane(self, left: QtWidgets.QWidget, right: QtWidgets.QWidget) -> QtWidgets.QSplitter:
            split = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
            split.addWidget(left); split.addWidget(right)
            split.setStretchFactor(0, 3); split.setStretchFactor(1, 2)
            return split

    # ------------------ Tab 1: Particles (lpp) ------------------
    class TabParticles(BaseTab):
        def __init__(self):
            super().__init__(bg_color="#f8fbff")
            # --- Top controls (with icons)
            self.btnBrowse=QtWidgets.QToolButton()
            self.btnBrowse.setIcon(load_icon("input.png", self.style().standardIcon(QtWidgets.QStyle.SP_DirOpenIcon)))
            self.btnClear=QtWidgets.QToolButton()
            self.btnClear.setIcon(load_icon("clear.png", self.style().standardIcon(QtWidgets.QStyle.SP_DialogResetButton)))
            self.outdirEdit=QtWidgets.QLineEdit(); self.outdirEdit.setPlaceholderText("Output folder (empty = same as input)")
            self.outdirBtn=QtWidgets.QToolButton()
            self.outdirBtn.setIcon(load_icon("output.png", self.style().standardIcon(QtWidgets.QStyle.SP_DirIcon)))
            hbar, hlabel = self.build_header_bar("Particles Dump → VTK", [self.btnBrowse, self.btnClear, self.outdirEdit, self.outdirBtn])
            self._header_label = hlabel

            # 詳細設定（折り畳み）
            self.det=Collapsible("Advanced")
            dlay=QtWidgets.QFormLayout(); dlay.setContentsMargins(4,4,4,4); dlay.setHorizontalSpacing(12); dlay.setVerticalSpacing(6)
            self.formatBox=QtWidgets.QComboBox(); self.formatBox.addItems(["ascii","binary"])
            self.backendBox=QtWidgets.QComboBox(); self.backendBox.addItems(["legacy","vtk"])
            self.cpuSpin=QtWidgets.QSpinBox(); self.cpuSpin.setRange(1, max(1, mp.cpu_count()))
            self.chunkSpin=QtWidgets.QSpinBox(); self.chunkSpin.setRange(1, 64)
            self.noOverwrite=QtWidgets.QCheckBox("--no-overwrite (skip existing .vtk)")
            # ラベル保持（i18n のため）
            self.lblFormat = QtWidgets.QLabel("format")
            self.lblBackend = QtWidgets.QLabel("backend")
            self.lblCPUs = QtWidgets.QLabel("CPUs")
            self.lblChunk = QtWidgets.QLabel("chunksize")
            dlay.addRow(self.lblFormat, self.formatBox)
            dlay.addRow(self.lblBackend, self.backendBox)
            dlay.addRow(self.lblCPUs, self.cpuSpin)
            dlay.addRow(self.lblChunk, self.chunkSpin)
            dlay.addRow("", self.noOverwrite)
            self.det.setContentLayout(dlay)
            # 詳細設定ボタンツールチップ（タブ固有の説明）
            self.det.setToolTip("粒子 VTK 出力に関する詳細設定です。書式（ASCII/BINARY）、バックエンド、並列数、チャンクサイズ、既存ファイルのスキップを切替できます。")

            # Left (inputs)
            left = QtWidgets.QWidget()
            llay = QtWidgets.QVBoxLayout(left); llay.setContentsMargins(6,6,6,6); llay.setSpacing(8)
            self.files.setPlaceholder("Drag & Drop here (particle dump)")
            llay.addWidget(self.files, 1)
            llay.addWidget(self.det, 0)

            # Right (progress + log + run)
            right = QtWidgets.QWidget()
            rlay = QtWidgets.QVBoxLayout(right); rlay.setContentsMargins(6,6,6,6); rlay.setSpacing(8)
            self.runBtn=QtWidgets.QPushButton("Export VTK ▶"); self.runBtn.setProperty("accent","true"); self.runBtn.setCursor(QtCore.Qt.PointingHandCursor)
            self.runBtn.setIcon(load_icon("convert.png"))
            self.runBtn.setIconSize(QtCore.QSize(48,48))  # convert.png は大きく
            rlay.addWidget(self.progress, 0)
            rlay.addWidget(self.log, 1)
            rlay.addWidget(self.runBtn, 0, alignment=QtCore.Qt.AlignRight)

            # Main layout (header + splitter)
            lay=QtWidgets.QVBoxLayout(self); lay.setContentsMargins(8,8,8,8); lay.setSpacing(8)
            lay.addWidget(hbar, 0)
            lay.addWidget(self.two_pane(left, right), 1)

            # Connect
            self.btnBrowse.clicked.connect(self._browse)
            self.btnClear.clicked.connect(self.clear_files)
            self.outdirBtn.clicked.connect(self._pick_outdir)
            self.runBtn.clicked.connect(self._run)
            self.requestOutdirSet.connect(self.outdirEdit.setText)
            # ファイル追加のたびに自動推奨へ合わせる
            self.files.filesDropped.connect(lambda _p: self._retune_parallel_defaults())

            # 初期の自動推奨
            self._retune_parallel_defaults()

            # Apply initial i18n
            self.apply_i18n()
            # React to language change
            i18n.langChanged.connect(lambda _l: self.apply_i18n())

        def _retune_parallel_defaults(self):
            cp = guess_workers()
            ch = guess_chunksize(self.files.count(), cp)
            self.cpuSpin.setValue(cp)
            self.chunkSpin.setValue(ch)
            # 動的ツールチップ
            _ = i18n.tr
            self.cpuSpin.setToolTip(_("Recommended from your CPUs: {}", "お使いのCPUからの推奨値: {}").format(cp))
            self.chunkSpin.setToolTip(_("Recommended from files and CPUs: {}", "ファイル数とCPUからの推奨値: {}").format(ch))

        def apply_i18n(self):
            _ = i18n.tr
            self._header_label.setText(_("Particles Dump → VTK", "粒子 dump → VTK"))
            self.btnBrowse.setToolTip(_("Select dump files (pick one to collect numbered siblings)", "dumpファイルを選択（1つ選ぶと連番をまとめて取り込みます）"))
            self.btnClear.setToolTip(_("Clear list", "リストをクリア"))
            self.outdirEdit.setPlaceholderText(_("Output folder (empty = same as input)", "出力フォルダ（空なら入力と同じ）"))
            self.outdirBtn.setToolTip(_("Choose output folder", "出力フォルダを選択"))
            self.det.setTitle(_("Advanced settings", "詳細設定"))
            self.noOverwrite.setText(_("--no-overwrite (skip existing .vtk)", "--no-overwrite（既存 .vtk をスキップ）"))
            self.files.setPlaceholder(_("Drag & Drop here (particle dump)", "ここにドラッグ＆ドロップ（粒子 dump）"))
            self.runBtn.setText(_("Export VTK ▶", "VTKを書き出す ▶"))
            self.progress.setFormat(_("%p% in progress", "%p% 実施中"))
            # ラベル翻訳
            self.lblFormat.setText(_("format", "書式（VTKレガシー）"))
            self.lblBackend.setText(_("backend", "バックエンド"))
            self.lblCPUs.setText(_("CPUs", "並列数"))
            self.lblChunk.setText(_("chunksize", "チャンクサイズ"))
            # tooltips（詳細説明）
            self.formatBox.setToolTip(_("VTK legacy file encoding.\nASCII = human-readable, BINARY = smaller/faster.",
                                        "VTKレガシー出力のエンコーディング。\nASCII=読みやすい / BINARY=小さく高速"))
            self.backendBox.setToolTip(_("Writer backend.\n'legacy' = pure Python, 'vtk' = python-vtk required.",
                                         "書き出しバックエンド。\n'legacy'=純Python / 'vtk'=python-vtk が必要"))
            self.noOverwrite.setToolTip(_("Skip when target .vtk already exists.", "出力先の .vtk が存在する場合はスキップします。"))
            # タブ固有の詳細説明（前段で英語/日本語を直接設定済み）
            # 明示的にトグルに対しても再設定
            self.det.setToolTip(_("Advanced settings for particle conversion (format/backend/parallel/chunk/skip).",
                                  "粒子変換の詳細設定（書式/バックエンド/並列/チャンク/スキップ）"))

        def _browse(self):
            _ = i18n.tr
            files,_sel = QtWidgets.QFileDialog.getOpenFileNames(self, _("Select dump files", "dumpファイルを選択"), "", _("All Files (*)", "すべてのファイル (*)"))
            if not files: return
            expanded = expand_series(files)
            self.add_files(expanded)
            self.outdirEdit.setText(ensure_outdir_default_from_inputs(expanded))
            self._retune_parallel_defaults()
        def _pick_outdir(self):
            _ = i18n.tr
            d=QtWidgets.QFileDialog.getExistingDirectory(self, _("Choose output folder", "出力フォルダを選択"), self.outdirEdit.text() or os.getcwd())
            if d: self.outdirEdit.setText(d if d.endswith(os.sep) else d+os.sep)
        def _run(self):
            _ = i18n.tr
            files=self.take_all_files()
            if not files:
                self.gui_log(_("No inputs.\n", "入力がありません。\n"), "#ef4444"); return
            fmt=self.formatBox.currentText(); backend=self.backendBox.currentText()
            cpunum=int(self.cpuSpin.value()); chunksize=int(self.chunkSpin.value()); skip=self.noOverwrite.isChecked()
            outroot = self.outdirEdit.text() or ensure_outdir_default_from_inputs(files)
            self.log.clear(); self.progress.setValue(0)
            def fn(filelist, output_root, fmt, backend, cpunum, chunksize, skip, log, progress):
                log(_("Started...\n", "開始しました…\n"))
                run_lpp_parallel(filelist, output_root, fmt=fmt, backend=backend, cpunum=cpunum, chunksize=chunksize, no_overwrite=skip, log=log, progress=progress)
            self.runner=Runner(fn, (files, outroot, fmt, backend, cpunum, chunksize, skip), {})
            self.runner.logged.connect(self.gui_log); self.runner.progressed.connect(self.set_progress)
            self.runner.finished_ok.connect(lambda: self.gui_log(_("Done\n", "完了\n"), "#16a34a"))
            self.runner.failed.connect(lambda m: self.gui_log(_("Error: ", "エラー: ") + m + "\n", "#ef4444"))
            self.runner.start()

    # ------------------ Tab 2: Header rename ------------------
    class TabHeader(BaseTab):
        DEFAULT_HEADER="ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz"
        def __init__(self):
            super().__init__(bg_color="#f7fff7")
            # header bar
            self.btnBrowse=QtWidgets.QToolButton(); 
            self.btnBrowse.setIcon(load_icon("input.png", self.style().standardIcon(QtWidgets.QStyle.SP_DirOpenIcon)))
            self.btnClear=QtWidgets.QToolButton(); 
            self.btnClear.setIcon(load_icon("clear.png", self.style().standardIcon(QtWidgets.QStyle.SP_DialogResetButton)))
            hbar, hlabel = self.build_header_bar("Header convert ('ITEM: ENTRIES ...' replacement)", [self.btnBrowse, self.btnClear])
            self._header_label = hlabel

            # 編集部
            self.formCard = QtWidgets.QGroupBox("Replacement header")
            self.headerEdit=QtWidgets.QLineEdit(self.DEFAULT_HEADER)
            fl = QtWidgets.QFormLayout(self.formCard); fl.setContentsMargins(10,6,10,8)
            self.lblEntries = QtWidgets.QLabel("ENTRIES")
            fl.addRow(self.lblEntries, self.headerEdit)

            # 詳細設定
            self.det=Collapsible("Advanced settings")
            dlay=QtWidgets.QFormLayout(); dlay.setContentsMargins(4,4,4,4); dlay.setHorizontalSpacing(12); dlay.setVerticalSpacing(6)
            self.inplace=QtWidgets.QCheckBox("Overwrite in-place"); self.inplace.setChecked(False)
            self.suffixEdit=QtWidgets.QLineEdit("-renamed"); self.encEdit=QtWidgets.QLineEdit("utf-8")
            self.lblSuffix = QtWidgets.QLabel("suffix")
            self.lblEncoding = QtWidgets.QLabel("encoding")
            dlay.addRow(self.lblSuffix, self.suffixEdit)
            dlay.addRow(self.lblEncoding, self.encEdit)
            dlay.addRow("", self.inplace)
            self.det.setContentLayout(dlay)
            # タブ固有の説明
            self.det.setToolTip("ヘッダー変換（ITEM: ENTRIES ... の置換）の詳細設定です。上書き保存/サフィックス/文字コードを切替できます。")

            # Left
            left = QtWidgets.QWidget()
            llay = QtWidgets.QVBoxLayout(left); llay.setContentsMargins(6,6,6,6); llay.setSpacing(8)
            self.files.setPlaceholder("Drag & Drop here (text dump)")
            llay.addWidget(self.files,1)
            llay.addWidget(self.formCard,0)
            llay.addWidget(self.det,0)

            # Right
            right = QtWidgets.QWidget()
            rlay = QtWidgets.QVBoxLayout(right); rlay.setContentsMargins(6,6,6,6); rlay.setSpacing(8)
            self.runBtn=QtWidgets.QPushButton("Apply header replacement ▶"); self.runBtn.setProperty("accent","true")
            self.runBtn.setIcon(load_icon("convert.png"))
            self.runBtn.setIconSize(QtCore.QSize(48,48))
            rlay.addWidget(self.progress, 0)
            rlay.addWidget(self.log, 1)
            rlay.addWidget(self.runBtn, 0, alignment=QtCore.Qt.AlignRight)

            # Main
            lay=QtWidgets.QVBoxLayout(self); lay.setContentsMargins(8,8,8,8); lay.setSpacing(8)
            lay.addWidget(hbar, 0)
            lay.addWidget(self.two_pane(left, right), 1)

            # Connect
            self.btnBrowse.clicked.connect(self._browse)
            self.btnClear.clicked.connect(self.clear_files)
            self.runBtn.clicked.connect(self._run)

            # i18n
            self.apply_i18n()
            i18n.langChanged.connect(lambda _l: self.apply_i18n())

        def apply_i18n(self):
            _ = i18n.tr
            self._header_label.setText(_("Header convert ('ITEM: ENTRIES ...' replacement)", "ヘッダー変換（'ITEM: ENTRIES ...' 置換）"))
            self.btnBrowse.setToolTip(_("Open target files", "対象ファイルを開く"))
            self.btnClear.setToolTip(_("Clear list", "リストをクリア"))
            self.formCard.setTitle(_("Replacement header", "置換ヘッダー"))
            self.lblEntries.setText(_("ENTRIES", "ENTRIES"))
            self.det.setTitle(_("Advanced settings", "詳細設定"))
            self.lblSuffix.setText(_("suffix", "サフィックス"))
            self.lblEncoding.setText(_("encoding", "文字コード"))
            self.files.setPlaceholder(_("Drag & Drop here (text dump)", "ここにドラッグ＆ドロップ（テキストダンプ）"))
            self.runBtn.setText(_("Apply header replacement ▶", "ヘッダー置換を実行 ▶"))
            self.progress.setFormat(_("%p% in progress", "%p% 実施中"))
            # tooltips
            self.headerEdit.setToolTip(_("Header line that replaces 'ITEM: ENTRIES ...'", "置換する 'ITEM: ENTRIES ...' の行を指定します。"))
            self.inplace.setToolTip(_("Overwrite the original file. If unchecked, writes a new file with suffix.", "原本を上書きします。未チェックの場合はサフィックス付きで新規保存します。"))
            self.suffixEdit.setToolTip(_("Suffix of new file when not overwriting.", "上書きしない場合の出力ファイル名サフィックス。"))
            self.encEdit.setToolTip(_("Text encoding to read/write.", "入出力に用いるテキストエンコーディング。"))
            self.det.setToolTip(_("Advanced settings for header replacement (in-place/suffix/encoding).",
                                   "ヘッダー置換の詳細設定（上書き/サフィックス/文字コード）。"))

        def _browse(self):
            _ = i18n.tr
            files,_sel=QtWidgets.QFileDialog.getOpenFileNames(self, _("Select files", "ファイルを選択"), "", _("All Files (*)", "すべてのファイル (*)"))
            if not files: return
            expanded = expand_series(files)
            self.add_files(expanded)
        def _run(self):
            _ = i18n.tr
            files=self.take_all_files()
            if not files: self.gui_log(_("No inputs.\n", "入力がありません。\n"), "#ef4444"); return
            header=self.headerEdit.text().strip()
            inplace=self.inplace.isChecked(); suffix=self.suffixEdit.text(); encoding=self.encEdit.text() or "utf-8"
            self.log.clear(); self.progress.setValue(0)
            def fn(filelist, header, inplace, suffix, encoding, log, progress):
                n=len(filelist); log(_("Started...\n", "開始しました…\n"))
                for i,f in enumerate(sorted(set(filelist))):
                    out = replace_header(Path(f), header, inplace, suffix, encoding)
                    log(f"✅ {f} -> {out}\n", "#22c55e")
                    progress(i+1, n)
            self.runner=Runner(fn, (files, header, inplace, suffix, encoding), {})
            self.runner.logged.connect(self.gui_log); self.runner.progressed.connect(self.set_progress)
            self.runner.finished_ok.connect(lambda: self.gui_log(_("Done\n", "完了\n"), "#16a34a"))
            self.runner.failed.connect(lambda m: self.gui_log(_("Error: ", "エラー: ") + m + "\n", "#ef4444"))
            self.runner.start()

    # ------------------ Tab 3: Force network ------------------
    class TabForce(BaseTab):
        def __init__(self):
            super().__init__(bg_color="#fbf7ff")
            # header
            self.btnBrowse=QtWidgets.QToolButton()
            self.btnBrowse.setIcon(load_icon("input.png", self.style().standardIcon(QtWidgets.QStyle.SP_DirOpenIcon)))
            self.btnClear=QtWidgets.QToolButton()
            self.btnClear.setIcon(load_icon("clear.png", self.style().standardIcon(QtWidgets.QStyle.SP_DialogResetButton)))
            self.outdirEdit=QtWidgets.QLineEdit(); self.outdirEdit.setPlaceholderText("Output folder (empty = same as input)")
            self.outdirBtn=QtWidgets.QToolButton()
            self.outdirBtn.setIcon(load_icon("output.png", self.style().standardIcon(QtWidgets.QStyle.SP_DirIcon)))
            hbar, hlabel = self.build_header_bar("Force chain → VTK/VTU + Louvain", [self.btnBrowse, self.btnClear, self.outdirEdit, self.outdirBtn])
            self._header_label = hlabel

            # 詳細設定（すべてここに格納）
            self.det=Collapsible("Advanced settings")
            dlay=QtWidgets.QFormLayout(); dlay.setContentsMargins(4,4,4,4); dlay.setHorizontalSpacing(12); dlay.setVerticalSpacing(6)
            self.vtkFormat=QtWidgets.QComboBox(); self.vtkFormat.addItems(["vtk","vtu"])  # 既定: vtk
            self.encoding=QtWidgets.QComboBox(); self.encoding.addItems(["ascii","binary"])
            self.keepPeriodic=QtWidgets.QCheckBox("Keep periodic contacts")
            self.resolution=QtWidgets.QDoubleSpinBox(); self.resolution.setDecimals(3); self.resolution.setRange(0.001, 100.0); self.resolution.setValue(1.0)
            self.seed=QtWidgets.QSpinBox(); self.seed.setRange(0, 2**31-1); self.seed.setValue(42)
            self.pointdata=QtWidgets.QCheckBox("Also write PointData (node_community/degree/force_sum)")
            self.nanfill=QtWidgets.QDoubleSpinBox(); self.nanfill.setDecimals(3); self.nanfill.setRange(-1e9,1e9); self.nanfill.setValue(0.0)
            self.cpuSpin=QtWidgets.QSpinBox(); self.cpuSpin.setRange(1, max(1, mp.cpu_count()))
            self.chunkSpin=QtWidgets.QSpinBox(); self.chunkSpin.setRange(1, 64)
            self.noOverwrite=QtWidgets.QCheckBox("--no-overwrite (skip existing)")
            # ラベル保持（i18n）
            self.lblFmt = QtWidgets.QLabel("vtk-format")
            self.lblEnc = QtWidgets.QLabel("encoding")
            self.lblKeep = QtWidgets.QLabel("keep-periodic")
            self.lblRes = QtWidgets.QLabel("resolution")
            self.lblSeed = QtWidgets.QLabel("seed")
            self.lblPD = QtWidgets.QLabel("write-pointdata")
            self.lblNan = QtWidgets.QLabel("nan-fill (VTK legacy)")
            self.lblCPU = QtWidgets.QLabel("CPUs")
            self.lblChunk = QtWidgets.QLabel("chunksize")
            dlay.addRow(self.lblFmt, self.vtkFormat)
            dlay.addRow(self.lblEnc, self.encoding)
            dlay.addRow(self.lblKeep, self.keepPeriodic)
            dlay.addRow(self.lblRes, self.resolution)
            dlay.addRow(self.lblSeed, self.seed)
            dlay.addRow(self.lblPD, self.pointdata)
            dlay.addRow(self.lblNan, self.nanfill)
            dlay.addRow(self.lblCPU, self.cpuSpin)
            dlay.addRow(self.lblChunk, self.chunkSpin)
            dlay.addRow("", self.noOverwrite)
            self.det.setContentLayout(dlay)
            # タブ固有の説明
            self.det.setToolTip("force chain 変換の詳細設定です。出力形式、エンコーディング、周期接触、Louvain の分解能/シード、PointData 出力、NaN置換、並列・チャンク、既存スキップを切替できます。")

            # Left
            left = QtWidgets.QWidget()
            llay = QtWidgets.QVBoxLayout(left); llay.setContentsMargins(6,6,6,6); llay.setSpacing(8)
            self.files.setPlaceholder("Drag & Drop here (force chain dump)")
            llay.addWidget(self.files,1)
            # 操作説明
            self.hint = QtWidgets.QLabel("Selecting one file auto-collects numbered siblings (same as particles tab)"); self.hint.setProperty("hint","true")
            llay.addWidget(self.hint, 0)
            llay.addWidget(self.det,0)

            # Right
            right = QtWidgets.QWidget()
            rlay = QtWidgets.QVBoxLayout(right); rlay.setContentsMargins(6,6,6,6); rlay.setSpacing(8)
            self.runBtn=QtWidgets.QPushButton("Export network ▶"); self.runBtn.setProperty("accent","true")
            self.runBtn.setIcon(load_icon("convert.png"))
            self.runBtn.setIconSize(QtCore.QSize(48,48))
            rlay.addWidget(self.progress, 0)
            rlay.addWidget(self.log, 1)
            rlay.addWidget(self.runBtn, 0, alignment=QtCore.Qt.AlignRight)

            # Main
            lay=QtWidgets.QVBoxLayout(self); lay.setContentsMargins(8,8,8,8); lay.setSpacing(8)
            lay.addWidget(hbar, 0)
            lay.addWidget(self.two_pane(left, right), 1)

            # Connect
            self.btnBrowse.clicked.connect(self._browse)
            self.btnClear.clicked.connect(self.clear_files)
            self.outdirBtn.clicked.connect(self._pick_outdir)
            self.runBtn.clicked.connect(self._run)
            self.requestOutdirSet.connect(self.outdirEdit.setText)
            self.files.filesDropped.connect(lambda _p: self._retune_parallel_defaults())

            # 初期の自動推奨
            self._retune_parallel_defaults()

            # i18n
            self.apply_i18n()
            i18n.langChanged.connect(lambda _l: self.apply_i18n())

        def _retune_parallel_defaults(self):
            cp = guess_workers()
            ch = guess_chunksize(self.files.count(), cp)
            self.cpuSpin.setValue(cp)
            self.chunkSpin.setValue(ch)
            _ = i18n.tr
            self.cpuSpin.setToolTip(_("Recommended from your CPUs: {}", "お使いのCPUからの推奨値: {}").format(cp))
            self.chunkSpin.setToolTip(_("Recommended from files and CPUs: {}", "ファイル数とCPUからの推奨値: {}").format(ch))

        def apply_i18n(self):
            _ = i18n.tr
            self._header_label.setText(_("Force chain → VTK/VTU + Louvain", "force chain → VTK/VTU + Louvain"))
            self.btnBrowse.setToolTip(_("Select pair/local dump", "pair/local dump を選択"))
            self.btnClear.setToolTip(_("Clear list", "リストをクリア"))
            self.outdirEdit.setPlaceholderText(_("Output folder (empty = same as input)", "出力フォルダ（空なら入力と同じ）"))
            self.outdirBtn.setToolTip(_("Choose output folder", "出力フォルダを選択"))
            self.det.setTitle(_("Advanced settings", "詳細設定"))
            self.keepPeriodic.setText(_("Keep periodic contacts", "周期接触を残す"))
            self.pointdata.setText(_("Also write PointData (node_community/degree/force_sum)","PointDataも書き出す（node_community/degree/force_sum）"))
            self.noOverwrite.setText(_("--no-overwrite (skip existing)", "--no-overwrite（既存をスキップ）"))
            self.files.setPlaceholder(_("Drag & Drop here (force chain dump)", "ここにドラッグ＆ドロップ（force chain dump）"))
            self.hint.setText(_("Selecting one file auto-collects numbered siblings (same as particles tab)","1ファイル選択で同名プレフィクスの連番を自動収集（粒子タブと同様）"))
            self.runBtn.setText(_("Export network ▶", "ネットワークを書き出す ▶"))
            self.progress.setFormat(_("%p% in progress", "%p% 実施中"))
            # ラベル翻訳
            self.lblFmt.setText(_("vtk-format", "出力形式"))
            self.lblEnc.setText(_("encoding", "エンコーディング"))
            self.lblKeep.setText(_("keep-periodic", "周期接触を保持"))
            self.lblRes.setText(_("resolution", "分解能（Louvain）"))
            self.lblSeed.setText(_("seed", "乱数シード"))
            self.lblPD.setText(_("write-pointdata", "PointDataを書き出す"))
            self.lblNan.setText(_("nan-fill (VTK legacy)", "NaN/Inf 置換値（VTKレガシー）"))
            self.lblCPU.setText(_("CPUs", "並列数"))
            self.lblChunk.setText(_("chunksize", "チャンクサイズ"))
            # tooltips
            self.vtkFormat.setToolTip(_("Output format: 'vtk' (legacy polydata) or 'vtu' (XML).", "出力形式：'vtk'(レガシーPOLYDATA) または 'vtu'(XML)。"))
            self.encoding.setToolTip(_("ASCII is human-readable; BINARY is compact & fast.", "ASCII は読みやすく、BINARY はサイズが小さく高速です。"))
            self.keepPeriodic.setToolTip(_("If checked, keep periodic contacts as edges.", "チェックすると周期境界に跨る接触を残します。"))
            self.resolution.setToolTip(_("Resolution parameter for Louvain; higher -> more communities.", "Louvain の分解能。大きいほどコミュニティ数が増えやすいです。"))
            self.seed.setToolTip(_("Random seed for Louvain community detection.", "Louvain の乱数シード。"))
            self.pointdata.setToolTip(_("Also write node_community/degree/force_sum as PointData.", "node_community/degree/force_sum を PointData として出力します。"))
            self.nanfill.setToolTip(_("Value used to replace NaN/Inf when writing legacy VTK.", "VTK レガシー書き出し時に NaN/Inf を置換する値。"))
            self.noOverwrite.setToolTip(_("Skip writing when destination file already exists.", "出力先が存在する場合はスキップします。"))
            self.det.setToolTip(_("Advanced settings for force-network conversion (format/encoding/periodic/Louvain/PointData/Nan/parallel/chunk/skip).",
                                  "force chain 変換の詳細設定（出力形式/エンコード/周期接触/Louvain/PointData/NaN/並列/チャンク/スキップ）。"))

        def _browse(self):
            _ = i18n.tr
            files,_sel=QtWidgets.QFileDialog.getOpenFileNames(self, _("Select pair/local dump", "pair/local dump を選択"), "", _("All Files (*)", "すべてのファイル (*)"))
            if not files: return
            expanded = expand_series(files)  # ★ 連番自動収集（数字が中間でもOK）
            self.add_files(expanded)
            self.outdirEdit.setText(ensure_outdir_default_from_inputs(expanded))
            self._retune_parallel_defaults()
        def _pick_outdir(self):
            _ = i18n.tr
            d=QtWidgets.QFileDialog.getExistingDirectory(self, _("Choose output folder", "出力フォルダを選択"), self.outdirEdit.text() or os.getcwd())
            if d: self.outdirEdit.setText(d)
        def _run(self):
            _ = i18n.tr
            files=self.take_all_files()
            if not files:
                self.gui_log(_("No inputs.\n", "入力がありません。\n"), "#ef4444"); return
            outdir=self.outdirEdit.text() or ensure_outdir_default_from_inputs(files)
            args=dict(vtk_format=self.vtkFormat.currentText(),
                      encoding=self.encoding.currentText(),
                      keep_periodic=self.keepPeriodic.isChecked(),
                      resolution=float(self.resolution.value()),
                      seed=int(self.seed.value()),
                      write_pointdata=self.pointdata.isChecked(),
                      outdir=outdir,
                      nan_fill=float(self.nanfill.value()),
                      skip_existing=self.noOverwrite.isChecked(),
                      cpunum=int(self.cpuSpin.value()),
                      chunksize=int(self.chunkSpin.value()))
            self.log.clear(); self.progress.setValue(0)
            self.runner=Runner(run_force_batch, (files,), args)
            self.runner.logged.connect(self.gui_log); self.runner.progressed.connect(self.set_progress)
            self.runner.finished_ok.connect(lambda: self.gui_log(_("Done\n", "完了\n"), "#16a34a"))
            self.runner.failed.connect(lambda m: self.gui_log(_("Error: ", "エラー: ") + m + "\n", "#ef4444"))
            self.runner.start()

    # ------------------ Main Window ------------------
    class Main(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("D2V — dump2vtk: Convert LAMMPS and LIGGGHTS Dump Files to VTK")
            self.resize(1060, 680)
            self.setAcceptDrops(False)  # 入力リストのみドロップ受付

            # Central tabs
            self.tabs=QtWidgets.QTabWidget()
            self.tab1=TabParticles(); self.tab2=TabHeader(); self.tab3=TabForce()
            self.tabs.addTab(self.tab1, "Dump→VTK (Particles)")
            self.tabs.addTab(self.tab2, "Header Replace")
            self.tabs.addTab(self.tab3, "Dump→VTK (Force chain)")
            self.setCentralWidget(self.tabs)

            # Top toolbar
            tb = QtWidgets.QToolBar(); tb.setIconSize(QtCore.QSize(16,16)); tb.setMovable(False); tb.setFloatable(False)
            self.helpLabel = QtWidgets.QLabel("Drag & Drop files or use the folder icons")
            tb.addWidget(self.helpLabel)
            tb.addSeparator()
            # Language button with menu
            self.langBtn = QtWidgets.QToolButton()
            self.langBtn.setText("Language")
            self.langBtn.setPopupMode(QtWidgets.QToolButton.InstantPopup)
            self.langMenu = QtWidgets.QMenu(self.langBtn)
            self.langGroup = QtGui.QActionGroup(self.langMenu)
            self.actEn = QtGui.QAction("English", self.langMenu, checkable=True)
            self.actJa = QtGui.QAction("日本語", self.langMenu, checkable=True)
            self.langGroup.addAction(self.actEn); self.langGroup.addAction(self.actJa)
            self.actEn.setChecked(True)  # default English
            self.langMenu.addAction(self.actEn); self.langMenu.addAction(self.actJa)
            self.langBtn.setMenu(self.langMenu)
            tb.addWidget(self.langBtn)
            self.addToolBar(QtCore.Qt.TopToolBarArea, tb)

            # Connect language actions
            self.actEn.triggered.connect(lambda: i18n.set_lang("en"))
            self.actJa.triggered.connect(lambda: i18n.set_lang("ja"))
            i18n.langChanged.connect(self.apply_i18n)

            # Initial i18n
            self.apply_i18n()

        def apply_i18n(self, *_):
            _ = i18n.tr
            self.setWindowTitle(_("D2V — dump2vtk: Convert LAMMPS and LIGGGHTS Dump Files to VTK",
                                  "D2V — dump2vtk: LIGGGHTS/LAMMPS ダンプを VTK へ変換"))
            self.tabs.setTabText(0, _("Dump→VTK (Particles)", "Dump→VTK（粒子）"))
            self.tabs.setTabText(1, _("Header Replace", "ヘッダー変換"))
            self.tabs.setTabText(2, _("Dump→VTK (Force chain)", "Dump→VTK（force chain）"))
            self.helpLabel.setText(_("Drag & Drop files or use the folder icons",
                                     "ファイルをドラッグ＆ドロップ、または「フォルダ」アイコンから選択"))
            self.langBtn.setText(_("Language", "言語"))
            # Check menu state
            self.actEn.setChecked(i18n.lang=="en")
            self.actJa.setChecked(i18n.lang=="ja")

    # ------------------ Application bootstrap ------------------
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(GLOBAL_QSS)

    # --- Find D2V.ico and set as app/window icon ---
    if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
        base_dir = sys._MEIPASS
    else:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    icon_path = resource_path("D2V.ico")
    if os.path.exists(icon_path):
        icon = QtGui.QIcon(icon_path)
        app.setWindowIcon(icon)

    w = Main()
    if os.path.exists(icon_path):
        w.setWindowIcon(QtGui.QIcon(icon_path))
    w.show()
    return app.exec()

# =============================================================================
#                                   main
# =============================================================================

def main():
    parser=build_cli()
    args=parser.parse_args()
    if not args.cmd:
        # launch GUI
        return run_gui()
    else:
        return run_cli(args)

if __name__=="__main__":
    mp.freeze_support()
    sys.exit(main())
