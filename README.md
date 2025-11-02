# dump2vtk (English)

A single-file **GUI + CLI** tool (PySide6) to quickly convert LIGGGHTS/LAMMPS dumps to **VTK/VTU**.  
It supports **particle dumps**, **pair/local (force chain) dumps**, and a **header replacer** in one script.

- Highlights: multiprocessing & chunking, **streaming** read for force chains, auto-collection of numbered siblings (digits can appear in the **middle**), **Louvain** community detection and **community aggregates**.
- Minimal dependencies: `numpy` (GUI uses `PySide6`; optional: `vtk` / `networkx` / `python-louvain`).
- License: see [LICENCE](./LICENCE).

---

# Compiled version (windows)
The compiled version is distributed below.
<p align="center">
  <a href="https://opendemjapan.booth.pm/items/7611221">
    <img
      alt="Free download: dump2vtk_installer.zip (156 MB)"
      src="https://img.shields.io/badge/%E7%84%A1%E6%96%99%E3%83%80%E3%82%A6%E3%83%B3%E3%83%AD%E3%83%BC%E3%83%89-dump2vtk__installer.zip%20(156%20MB)-d32f2f?style=for-the-badge&logo=icloud&logoColor=white&labelColor=d32f2f">
  </a>
</p>


---

## 1. How to run with Python

### 1.1 Environment
- Python **3.8+** (recommended: 3.10 or newer)

### 1.2 Installation
```bash
# minimal
pip install numpy

# for the GUI
pip install PySide6

# optional features
pip install vtk networkx python-louvain
```

### 1.3 Start the GUI
```bash
python dump2vtk.py
```
- The GUI launches. Drag & drop files into the list and configure **Advanced** options as needed.

### 1.4 Run from the CLI
Three subcommands are available. Use `python dump2vtk.py ...` on all OSes.

**(a) Particles dump → VTK (legacy / python-vtk)**  
```bash
python dump2vtk.py lpp DUMP_FILES... \
  -o OUTROOT --format {ascii|binary} \
  --cpunum N --chunksize K --no-overwrite
```

**(b) Header replacement (`ITEM: ENTRIES ...`)**  
```bash
python dump2vtk.py rename \
  -H "ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz Fnx Fny Fnz Ftx Fty Ftz" \
  inputs*.dump --inplace
```

**(c) Force chain → VTK/VTU (Louvain + aggregates)**  
```bash
python dump2vtk.py force forcechain-*.dump \
  --vtk-format {vtk|vtu} --encoding {ascii|binary} --keep-periodic \
  --resolution 1.0 --seed 42 --write-pointdata \
  --outdir out/ --nan-fill 0.0 \
  --cpunum N --chunksize K --no-overwrite
```

---

## 2. GUI usage

The GUI has **three tabs** (see screenshots below).

### 2.1 Particles: Dump → VTK

![Particles tab — auto-collect series & parallel conversion](1.png)

1. Drag & drop **one** dump file into the list on the left.  
   The app **auto-collects numbered siblings** sharing the same prefix (digits can be at the end or **in the middle**).
2. Optionally set the **output folder** at the top-right (blank = same as input).
3. In **Advanced settings**, choose **ASCII/BINARY**, CPU count, chunk size, and `--no-overwrite` (skip existing outputs).  
   *Note*: You can switch the writer backend between **legacy (pure Python)** and **vtk (python-vtk bindings)**.
4. Click **Export VTK**. If the log shows  
   `Start: particle VTK (parallel)... Done.` the conversion succeeded.

---

### 2.2 Header Replace: replace `ITEM: ENTRIES ...`

![Header Replace tab — bulk replacement of ITEM: ENTRIES ...](2.png)

1. Drop pair/local-style dumps (e.g., `fc0.dump, fc2000.dump, ...`).  
   WSL paths like `\\wsl.localhost\...` are supported.
2. Enter the following one-line header in **Replacement header**:  
   ```text
   ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz Fnx Fny Fnz Ftx Fty Ftz
   ```
3. Click **Apply header replacement**. Each file is duplicated as `*-renamed.dump`.

---

### 2.3 Force chain: force network (VTK/VTU output)

![Force chain tab — batch processing with Louvain aggregates](3.png)

1. Select **one** `*-renamed.dump`; the app **auto-collects** its numbered siblings.
2. Click **Export network**. A **VTK/VTU** file is written per timestep.  
   When the log shows `Start: force chain batch (parallel)... All done.` the batch has finished.

---

## 3. VTK examples (ParaView)

Examples of the output files in ParaView.

![VTK example — particles](particle.png)  
*Particles (VTK POLYDATA)*

![VTK example — force chain](force_chain.png)  
*Force network (polyline edges)*

![VTK example — Louvain network](louvain.png)  
*Louvain-colored force network*

---

## 4. I/O essentials

- **Particle dumps**: detects `x/y/z` vs `xs/ys/zs` and **unscales** if needed. Consecutive `...x/..y/...z` and patterns like `f_*/c_*/v_*[1..3]` are treated as **vectors**; other columns become **scalars**.  
  The legacy VTK writer outputs `POINTS`, `VERTICES`, and `POINT_DATA`, plus a per-snapshot **bounding box** (RECTILINEAR_GRID). The python-vtk backend produces equivalent geometry via VTK bindings.
- **Force chain**: expects **12 required columns** in `ENTRIES`:  
  `x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz`  
  Additional columns are auto-detected as scalars/vectors. Edges come from `(id1, id2)`. Derived fields include `force = ||(fx, fy, fz)||` and `connectionLength`. **Louvain** adds `community` and `intra_comm`. With `--write-pointdata`, node-level `node_community`, `node_degree`, `node_force_sum` are included. For each attribute, **community mean/sum** is also emitted onto edges.

---

## 5. Auto-collection of numbered series

- **prefix + digits + suffix + extension** (supports `.gz`): adding a single file will discover and include its **sibling files**.  
- Works even when the digits are in the **middle** (e.g., `forcechain-16000-renamed.dmp`).  
- Wildcards like `forcechain-*.dmp` are also supported; duplicates are removed automatically.

---

## 6. Acknowledgments

- The design for dump handling and VTK output was inspired by **Pizza.py / LPP (LIGGGHTS post‑processing)**.  
- Examples are adapted from **compaction_LIGGGHTS**. Many thanks to the authors and community.

---

## 7. License

This repository includes code derived from **Pizza.py/LPP**, therefore the combined distribution is provided under **GPL‑2.0**.  
All new/original parts authored for this project are additionally offered under the **MIT License**. When distributing with GPL‑covered parts, **GPL‑2.0** terms take precedence.

- SPDX: `GPL-2.0-only AND MIT`  
- See [`LICENCE`](./LICENCE) for details.

---


