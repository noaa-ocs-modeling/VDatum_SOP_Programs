#!/usr/bin/env python3
"""
vpop_nc_Pacific_mp4.py
Vectorized Python port of vpop29nc4.f for multi-polygon VDatum processing.
Designed to run one polygon per SLURM array task.

Usage (single polygon test):
    export SLURM_ARRAY_TASK_ID=0
    python -u vpop_nc_Pacific_mp4.py

Usage (batch via SLURM):
    sbatch run_vpop_Pacific.sh

Environment variables:
    SLURM_ARRAY_TASK_ID   Task index set automatically by SLURM
    VPOP_INDIR            Directory containing vpop29_PA*_01.in files
                          Example default:
                          /work2/noaa/vdatum/.../vpop/out_vpop29_in_files
    MARINE_GRID_DIR       Directory containing marine grid NetCDF files
                          Example default:
                          /work2/noaa/vdatum/.../vgrid/out_marine_grid_nc
    VPOP_OUTDIR           Directory for polygon-specific output
    MESH_PATH             Shared ADCIRC mesh file (fort.14 or mesh NetCDF)
                          Required when datum inputs are NetCDF files that
                          contain values but not mesh geometry/connectivity
    VPOP_MODE             standard or svu
    VPOP_NWORKERS         Number of parallel workers
                          Default: SLURM_CPUS_PER_TASK or 4
    VPOP_K_QUERY          KDTree query size (optional, default: 512)
    VPOP_SMALL_EPS        Small-scale fill tolerance (optional, default: 1e-10)
    VPOP_SMALL_MAX_NODES  Max nodes for small-scale fill (optional, default: 0)

Input files:
    vpop29_PA*_01.in

Required Python packages:
    numpy scipy netCDF4 matplotlib

Optional:
    shapely (for polygon work)

Notes:
    - This is a practical, testable port rather than a line-by-line translation.
    - Expensive nested loops are replaced with NumPy, SciPy cKDTree,
      Matplotlib TriFinder, and Numba-accelerated routines.
    - When datum inputs are NetCDF files containing only node values,
      a shared ADCIRC mesh must be provided through MESH_PATH.
    - v3 mode: improved small-scale fill tolerance for mfill3/mfill4 parity.
    - icorr_ord (tidal-order correction) is disabled, matching vpop29nc4.f behavior.
    - icorr_dat (station correction) is not implemented.
    - MSL input is not supported (will raise error).

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional
import argparse
import math
import os
import sys
import glob
import numpy as np

try:
    from netCDF4 import Dataset
except Exception as exc:  # pragma: no cover
    raise SystemExit("Missing dependency: netCDF4. Install with conda/pip.") from exc

try:
    from scipy.spatial import cKDTree
except Exception as exc:  # pragma: no cover
    raise SystemExit("Missing dependency: scipy. Install scipy.") from exc

try:
    from numba import njit
    HAVE_NUMBA = True
except Exception:
    HAVE_NUMBA = False
    njit = None

try:
    import matplotlib.tri as mtri
    from matplotlib.path import Path as MplPath
except Exception as exc:  # pragma: no cover
    raise SystemExit("Missing dependency: matplotlib. Install matplotlib.") from exc

DEF_M8 = -88.8888
DEF_P9 = 9.999
RAD = 0.01745329


@dataclass
class VpopConfig:
    title: str
    icorr_ord: int
    icorr_dat: int
    ivoid_stop: int
    ismooth: int
    print_points: list[tuple[int, int]]
    ifldprt: list[int]
    value1: float
    value2: float
    fdatum_new: str
    msf_file: str
    nontidal: str
    fpolygon: str
    gpolygon: str
    fdatum_msl: str
    unfilled: str
    fitide: str
    vgridfile: str
    nfmax: int
    ngrids: int
    fdatum_in: list[list[str]]  # [nf][ng]
    fdatum_out: list[str]


@dataclass
class MarineGrid:
    isv: np.ndarray  # shape (nx, ny), 0 land, 1 water, 2+ layers
    x: np.ndarray    # shape (nx,)
    y: np.ndarray    # shape (ny,)
    dx: float
    dy: float
    box: np.ndarray


@dataclass
class AdcircData:
    x: np.ndarray             # nodes
    y: np.ndarray
    triangles: np.ndarray     # shape (nt, 3), zero based
    ds: np.ndarray            # shape (nodes, nfmax)
    iss: np.ndarray           # node inside BP/G polygon exactly
    iss2: np.ndarray          # expanded to triangles, excludes DEF_P9 nodes


@dataclass
class PolySegments:
    segments: list[np.ndarray]  # each array shape (n, 2), columns x,y
    poly_ids: list[int]

    @property
    def empty(self) -> bool:
        return len(self.segments) == 0


def _strip_line(line: str) -> str:
    return line.rstrip("\n").strip()


def _is_none(s: str) -> bool:
    return s.strip().lower() == "none"


def _resolve(path: str, base: Path) -> str:
    p = Path(path.strip())
    if p.is_absolute() or _is_none(path):
        return str(p) if not _is_none(path) else "none"
    return str((base / p).resolve())

def infer_vpop_mode(cfg: VpopConfig) -> str:
    joined = " ".join(cfg.fdatum_out).lower()
    if "_svu" in joined:
        return "svu"
    return "standard"

def wrap_lon_to_360(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x < 0.0, x + 360.0, x)

def wrap_lon_to_180(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x > 180.0, x - 360.0, x)
    
def parse_vpop_input(infile: str | Path) -> VpopConfig:
    """Parse the vpop29 input file written by the user's MATLAB prep script."""
    infile = Path(infile).resolve()
    base = infile.parent
    raw = infile.read_text(errors="ignore").splitlines()

    # Remove empty-only lines but keep the line order otherwise.
    lines = [l.rstrip("\n") for l in raw if l.strip()]
    k = 0

    def next_line() -> str:
        nonlocal k
        if k >= len(lines):
            raise ValueError(f"Unexpected end of input while parsing {infile}")
        s = lines[k]
        k += 1
        return s

    def first_tokens(line: str) -> list[str]:
        # vpop .in files commonly contain inline comments after "<", for example:
        #   0  72 292  46 65  1389 86  <nprt, i,jprnt
        # Strip those comments before tokenizing so we only parse numeric/path fields.
        return line.split("<", 1)[0].split()

    def first_field(line: str) -> str:
        toks = first_tokens(line)
        if not toks:
            raise ValueError(f"Expected a field before inline comment in line: {line!r}")
        return toks[0]

    title = next_line().strip()
    icorr_ord = int(first_tokens(next_line())[0])
    icorr_dat = int(first_tokens(next_line())[0])
    ivoid_stop = int(first_tokens(next_line())[0])
    ismooth = int(first_tokens(next_line())[0])

    toks = first_tokens(next_line())
    nprt = int(toks[0])
    pts = []
    vals = [int(v) for v in toks[1:]]
    for ii in range(min(nprt, len(vals) // 2)):
        pts.append((vals[2 * ii], vals[2 * ii + 1]))

    ifldprt = [int(v) for v in first_tokens(next_line())[:5]]

    toks = first_tokens(next_line())
    value1, value2, fdatum_new = float(toks[0]), float(toks[1]), toks[2]

    msf_file = first_field(next_line())
    nontidal = first_field(next_line())
    fpolygon = first_field(next_line())
    gpolygon = first_field(next_line())
    fdatum_msl = first_field(next_line())
    unfilled = first_field(next_line())
    fitide = first_field(next_line())
    vgridfile = first_field(next_line())

    toks = first_tokens(next_line())
    nfmax, ngrids = int(toks[0]), int(toks[1])

    fdatum_in: list[list[str]] = []
    fdatum_out: list[str] = []
    for _nf in range(nfmax):
        one = []
        for _ng in range(ngrids):
            one.append(first_field(next_line()))
        fdatum_in.append(one)
        fdatum_out.append(first_field(next_line()))

    return VpopConfig(
        title=title,
        icorr_ord=icorr_ord,
        icorr_dat=icorr_dat,
        ivoid_stop=ivoid_stop,
        ismooth=ismooth,
        print_points=pts,
        ifldprt=ifldprt,
        value1=value1,
        value2=value2,
        fdatum_new=_resolve(fdatum_new, base),
        msf_file=_resolve(msf_file, base),
        nontidal=_resolve(nontidal, base),
        fpolygon=_resolve(fpolygon, base),
        gpolygon=_resolve(gpolygon, base),
        fdatum_msl=_resolve(fdatum_msl, base),
        unfilled=_resolve(unfilled, base),
        fitide=_resolve(fitide, base),
        vgridfile=_resolve(vgridfile, base),
        nfmax=nfmax,
        ngrids=ngrids,
        fdatum_in=[[_resolve(p, base) for p in row] for row in fdatum_in],
        fdatum_out=[_resolve(p, base) for p in fdatum_out],
    )


def read_marine_grid_nc(filename: str | Path) -> MarineGrid:
    """Read vgridder/vpop marine grid NetCDF: box + field.

    IMPORTANT PARITY DETAIL
    -----------------------
    The NetCDF variable is stored on disk as shape (jmax, imax), which for
    AKai_de_01 is (3145, 3008).  Fortran works internally as field(i,j), i.e.
    shape (imax, jmax) = (3008, 3145).

    Therefore this reader transposes the NetCDF field immediately so the rest
    of the Python code uses the same internal order as Fortran:
        grid.isv.shape == (imax, jmax)
        grid.x.size    == imax
        grid.y.size    == jmax

    Writers transpose back to NetCDF order.
    """
    filename = Path(filename)
    with Dataset(filename, "r") as nc:
        box = np.asarray(nc.variables["box"][:], dtype=np.float64)
        field_nc = np.asarray(nc.variables["field"][:])

    if field_nc.ndim != 2:
        raise ValueError(f"Expected 2-D field in {filename}, got shape {field_nc.shape}")

    # NetCDF order is (jmax, imax); Fortran/Python internal order is (imax, jmax).
    field = field_nc.T.copy()

    imax, jmax = field.shape
    y0, x0, dy, dx = [float(v) for v in box]
    x = x0 + np.arange(imax, dtype=np.float64) * dx
    y = y0 + np.arange(jmax, dtype=np.float64) * dy
    return MarineGrid(isv=field.astype(np.int16), x=x, y=y, dx=dx, dy=dy, box=box)


def read_one_adcirc_ascii(filename: str | Path):
    """Read one ADCIRC-style triangular grid datum file.

    Format expected by vpop29nc4.f:
        title
        nt nodes
        node_id x y datum
        elem_id np n1 n2 n3
    Node IDs in the element section are 1-based in the source file.
    """
    filename = Path(filename)
    with filename.open("r", errors="ignore") as f:
        title = f.readline().rstrip("\n")
        nt, nodes = [int(v) for v in f.readline().split()[:2]]
        arr = np.loadtxt([next(f) for _ in range(nodes)], dtype=np.float64)
        elem = np.loadtxt([next(f) for _ in range(nt)], dtype=np.int64)

    if arr.ndim == 1:
        arr = arr[None, :]
    if elem.ndim == 1:
        elem = elem[None, :]
    x = arr[:, 1].astype(np.float64)
    y = arr[:, 2].astype(np.float64)
    val = arr[:, 3].astype(np.float64)
    np_col = elem[:, 1]
    if not np.all(np_col == 3):
        raise ValueError(f"Only 3-node ADCIRC elements are supported. {filename} has non-triangle elements.")
    tri = elem[:, 2:5].astype(np.int64) - 1
    return title, x, y, val, tri

def read_adcirc_full_mesh(adcirc_path: str):
    """Read shared ADCIRC fort.14 mesh geometry."""
    with open(adcirc_path, "rt") as f:
        _title = f.readline().rstrip("\n")
        parts = f.readline().split()
        nt, nodes = int(float(parts[0])), int(float(parts[1]))

        xs = np.zeros(nodes, dtype=np.float64)
        ys = np.zeros(nodes, dtype=np.float64)
        ds = np.zeros(nodes, dtype=np.float64)

        for _ in range(nodes):
            row = f.readline().split()
            nid = int(float(row[0])) - 1
            xs[nid] = float(row[1])
            ys[nid] = float(row[2])
            ds[nid] = float(row[3])

        tri = np.zeros((nt, 3), dtype=np.int32)
        for _ in range(nt):
            row = f.readline().split()
            eid = int(float(row[0])) - 1
            tri[eid, :] = [int(float(x)) - 1 for x in row[2:5]]

    return xs, ys, ds, tri


def read_adcirc_netcdf_mesh(mesh_path: str):
    """Read shared ADCIRC mesh geometry from NetCDF."""
    with Dataset(mesh_path, "r") as ncm:
        if "x" in ncm.variables:
            xs = np.asarray(ncm.variables["x"][:], dtype=np.float64)
        elif "lon" in ncm.variables:
            xs = np.asarray(ncm.variables["lon"][:], dtype=np.float64)
        else:
            raise KeyError(f"Mesh file {mesh_path} missing longitude variable (x/lon)")

        if "y" in ncm.variables:
            ys = np.asarray(ncm.variables["y"][:], dtype=np.float64)
        elif "lat" in ncm.variables:
            ys = np.asarray(ncm.variables["lat"][:], dtype=np.float64)
        else:
            raise KeyError(f"Mesh file {mesh_path} missing latitude variable (y/lat)")

        if "element" in ncm.variables:
            tri = np.asarray(ncm.variables["element"][:], dtype=np.int32) - 1
        elif "ele" in ncm.variables:
            tri = np.asarray(ncm.variables["ele"][:], dtype=np.int32) - 1
        else:
            raise KeyError(f"Mesh file {mesh_path} missing connectivity variable (element/ele)")

        if "depth" in ncm.variables:
            ds = np.asarray(ncm.variables["depth"][:], dtype=np.float64)
        elif "z" in ncm.variables:
            ds = np.asarray(ncm.variables["z"][:], dtype=np.float64)
        else:
            ds = np.zeros_like(xs)

    return np.array(xs).flatten(), np.array(ys).flatten(), np.array(ds).flatten(), tri


def read_pacific_value_nc(nc_path: str, mode: str = "standard") -> np.ndarray:
    """
    Read Pacific per-datum NetCDF values.

    standard mode -> d
    svu mode      -> diaPA_real if present, else diaPA real part, else zeros
    """
    with Dataset(nc_path, "r") as nc:
        if "d" not in nc.variables:
            raise KeyError(f"{nc_path} missing variable 'd'")
        d = np.asarray(nc.variables["d"][:], dtype=np.float64)

        if "diaPA_real" in nc.variables:
            diaPA = np.asarray(nc.variables["diaPA_real"][:], dtype=np.float64)
        elif "diaPA" in nc.variables:
            diaPA = nc.variables["diaPA"][:]
            if np.iscomplexobj(diaPA):
                diaPA = np.real(diaPA)
            diaPA = np.asarray(diaPA, dtype=np.float64)
        else:
            diaPA = np.zeros_like(d)

        if "diaPA_imag" in nc.variables:
            diaPA_imag = np.asarray(nc.variables["diaPA_imag"][:], dtype=np.float64)
            n_imag = int(np.count_nonzero(np.abs(diaPA_imag) > 0.0))
            if n_imag > 0:
                print(f"WARNING: {os.path.basename(nc_path)} has {n_imag} nodes with nonzero diaPA_imag", flush=True)

    mode = mode.lower()
    if mode in ("standard", "datum", "d"):
        return d
    elif mode in ("svu", "unc", "uncertainty", "diapa"):
        return diaPA
    else:
        raise ValueError(f"Unsupported VPOP mode: {mode}")
        
def read_adcirc_all_ascii(fdatum_in: list[list[str]], nfmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ASCII/GRD datum input path for meshes stored directly in the datum files."""
    all_x = []
    all_y = []
    all_tri = []
    ds_cols = []
    geom_nodes_total = None

    for nf in range(nfmax):
        vals_this = []
        node_offset = 0
        nf_x = []
        nf_y = []
        nf_tri = []
        for ng, path in enumerate(fdatum_in[nf]):
            _title, x, y, val, tri = read_one_adcirc_ascii(path)
            vals_this.append(val)
            if nf == 0:
                nf_x.append(x)
                nf_y.append(y)
                nf_tri.append(tri + node_offset)
            node_offset += len(x)

        vals_cat = np.concatenate(vals_this)
        ds_cols.append(vals_cat)

        if nf == 0:
            all_x = [np.concatenate(nf_x)]
            all_y = [np.concatenate(nf_y)]
            all_tri = [np.vstack(nf_tri)]
            geom_nodes_total = len(all_x[0])
        else:
            if len(vals_cat) != geom_nodes_total:
                raise ValueError(
                    f"Datum {nf+1} has {len(vals_cat)} nodes, expected {geom_nodes_total}. "
                    "vpop expects all datum files to share geometry."
                )

    x = all_x[0]
    y = all_y[0]
    tri = all_tri[0]
    ds = np.vstack(ds_cols).T.astype(np.float64)
    return x, y, tri, ds


def read_adcirc_all_pacific(
    fdatum_in: list[list[str]],
    nfmax: int,
    mesh_path: str,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pacific NetCDF path: shared mesh geometry + per-datum value NetCDF."""
    mesh_ext = Path(mesh_path).suffix.lower()
    if mesh_ext == ".nc":
        x, y, _mesh_ds, tri = read_adcirc_netcdf_mesh(mesh_path)
    else:
        x, y, _mesh_ds, tri = read_adcirc_full_mesh(mesh_path)

    geom_nodes_total = len(x)
    ds_cols = []

    for nf in range(nfmax):
        vals_this = []
        for ng, path in enumerate(fdatum_in[nf]):
            vals = read_pacific_value_nc(path, mode=mode)
            vals_this.append(vals)

        vals_cat = np.concatenate(vals_this)
        if len(vals_cat) != geom_nodes_total:
            raise ValueError(
                f"Datum {nf+1} has {len(vals_cat)} nodes, expected {geom_nodes_total}. "
                "Pacific vpop expects all datum files to match the shared mesh geometry."
            )
        ds_cols.append(vals_cat)

    ds = np.vstack(ds_cols).T.astype(np.float64)
    return x, y, tri, ds


def read_adcirc_all(
    fdatum_in: list[list[str]],
    nfmax: int,
    mode: str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Dispatcher:
      - .nc input  -> shared mesh + NetCDF values
      - otherwise  -> ASCII/GRD datum input path
    """
    first = str(fdatum_in[0][0]).lower()

    if first.endswith(".nc"):
        mesh_path = os.environ.get("MESH_PATH", "").strip()
        if not mesh_path:
            raise RuntimeError(
                "Pacific NetCDF datum input detected, but MESH_PATH is not set."
            )

        if mode is None or not mode.strip():
            mode = os.environ.get("VPOP_MODE", "standard").strip().lower() or "standard"

        print(f"  Pacific NetCDF mode: {mode}", flush=True)
        print(f"  Using shared mesh: {mesh_path}", flush=True)
        return read_adcirc_all_pacific(fdatum_in, nfmax, mesh_path, mode)

    return read_adcirc_all_ascii(fdatum_in, nfmax)
def detect_lon_lat(col0: np.ndarray, col1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Detect which column is longitude and which is latitude.

    Same idea as vgridder:
    - latitude-like values are mostly within [-90, 90]
    - longitude-like values may be outside that range, especially in 0..360
    """
    a = np.asarray(col0, dtype=np.float64)
    b = np.asarray(col1, dtype=np.float64)

    a_lat_like = np.mean(np.abs(a) <= 90.0)
    b_lat_like = np.mean(np.abs(b) <= 90.0)

    if a_lat_like > 0.95 and b_lat_like < 0.95:
        return b, a
    if b_lat_like > 0.95 and a_lat_like < 0.95:
        return a, b

    # fallback: keep original order
    return a, b    

def read_poly_segments(filename: str | Path | None) -> PolySegments:
    """Read BP/grid/non-tidal polygon DAT file in vpop rd_poly style.

    Updated for Pacific:
    - detect lon/lat columns like vgridder instead of relying only on a<0
    - still preserves segment start/end handling from original logic
    """
    if filename is None or _is_none(str(filename)):
        return PolySegments([], [])

    filename = Path(filename)
    if not filename.exists():
        raise FileNotFoundError(filename)

    rows = []
    with filename.open("r", errors="ignore") as f:
        for line in f:
            toks = line.split()
            if len(toks) < 3:
                continue
            try:
                nums = [float(t) for t in toks[:4]]
            except ValueError:
                continue
            rows.append(nums)

    if not rows:
        return PolySegments([], [])

    data = np.asarray(rows, dtype=np.float64)

    col0 = data[:, 0]
    col1 = data[:, 1]
    flags = np.rint(data[:, 2]).astype(np.int32)

    # Detect lon/lat like vgridder
    lon, lat = detect_lon_lat(col0, col1)

    # Polygon IDs: if 4th column exists and filename contains xyij, use it
    if ("xyij" in filename.name) and (data.shape[1] >= 4):
        pids_raw = np.rint(data[:, 3]).astype(np.int32)
    else:
        pids_raw = None

    # Build segment numbering like original code
    pids = []
    isum = 0
    for i in range(len(flags)):
        isum += int(flags[i])
        ns = (isum + 1) // 2
        pid = int(pids_raw[i]) if pids_raw is not None else ns
        pids.append(pid)

    starts = []
    isum = 0
    for idx, flag in enumerate(flags):
        isum += int(flag)
        if flag == 1 and (isum % 2) != 0:
            starts.append(idx)

    if not starts:
        starts = [0]
    starts.append(len(lon))

    segs = []
    seg_pids = []
    for s0, s1 in zip(starts[:-1], starts[1:]):
        pts = np.column_stack([lon[s0:s1], lat[s0:s1]])
        if len(pts) >= 3:
            segs.append(np.asarray(pts, dtype=np.float64))
            seg_pids.append(int(pids[s0]))

    return PolySegments(segs, seg_pids)

def contains_points(poly: PolySegments, x: np.ndarray, y: np.ndarray, first_segment_only: bool = False) -> np.ndarray:
    """Vectorized point-in-polygon using matplotlib.path.

    The Fortran has custom boundary perturbation. This uses Path.contains_points;
    a tiny positive radius treats boundary points as inside, which is closer to
    the vpop behavior for practical grid screening.
    """
    if poly.empty:
        return np.ones_like(x, dtype=bool)
    pts = np.column_stack([x.ravel(), y.ravel()])
    mask = np.zeros(pts.shape[0], dtype=bool)
    segments = poly.segments[:1] if first_segment_only else poly.segments
    for seg in segments:
        path = MplPath(seg, closed=True)
        mask |= path.contains_points(pts, radius=1.0e-12)
    return mask.reshape(np.shape(x))


def screen_adcirc_nodes(
    x: np.ndarray,
    y: np.ndarray,
    triangles: np.ndarray,
    ds: np.ndarray,
    bpoly: PolySegments,
    gpoly: PolySegments,
) -> tuple[np.ndarray, np.ndarray]:
    """Implement adcirc_inside_poly4 in vectorized form."""
    iss = np.ones(len(x), dtype=bool)
    if not bpoly.empty:
        iss &= contains_points(bpoly, x, y)
        if not np.any(iss):
            raise RuntimeError("No ADCIRC nodes are inside the bounding polygon.")
    if not gpoly.empty:
        iss &= contains_points(gpoly, x, y)
        if not np.any(iss):
            raise RuntimeError("No ADCIRC nodes are inside the grid polygon.")

    # If one node of a triangle is inside, all three triangle nodes become usable.
    tri_has_inside = np.any(iss[triangles], axis=1)
    iss2 = iss.copy()
    if np.any(tri_has_inside):
        iss2[np.unique(triangles[tri_has_inside].ravel())] = True

    # Fortran only checks ds(:,1) for 9.999 when disabling nodes.
    undefined = np.abs(ds[:, 0] - DEF_P9) <= 0.001
    iss2[undefined] = False
    return iss.astype(np.int8), iss2.astype(np.int8)


def initialize_fields(grid: MarineGrid, nfmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nx, ny = grid.isv.shape
    nfp = nfmax + 1
    dv = np.full((nx, ny, nfp), DEF_M8, dtype=np.float64)
    mfill = np.zeros((nx, ny), dtype=np.int16)
    itide = np.where(grid.isv >= 1, 1, 0).astype(np.int16)
    return dv, mfill, itide


def get_weights(grid: MarineGrid) -> np.ndarray:
    """Fortran get_weights."""
    nx, ny = grid.isv.shape
    mid_j = ny // 2
    mid_i = nx // 2
    ddx = math.cos(RAD * grid.y[mid_j])
    ddy = (grid.y[min(mid_j + 1, ny - 1)] - grid.y[mid_j]) / (grid.x[min(mid_i + 1, nx - 1)] - grid.x[mid_i])
    power = 1.0
    wt = np.empty(8, dtype=np.float64)
    wt[0:2] = 1.0 / (ddx ** power)       # left/right
    wt[2:4] = 1.0 / (ddy ** power)       # down/up
    wt[4:8] = 1.0 / (math.sqrt(ddx * ddx + ddy * ddy) ** power)  # diagonals
    return wt


def mark_non_tidal(
    grid: MarineGrid,
    bpoly: PolySegments,
    ntpoly: PolySegments,
    dv: np.ndarray,
    mfill: np.ndarray,
    itide: np.ndarray,
    value1: float,
    value2: float,
    nfmax: int,
) -> int:
    """Vectorized non_tidals/get_itide behavior."""
    if ntpoly.empty:
        # Fortran still initializes mfill for all water/layers later.
        return 0

    X, Y = np.meshgrid(grid.x, grid.y, indexing="ij")
    water = grid.isv >= 1
    in_bp_seg1 = contains_points(bpoly, X, Y, first_segment_only=True) if not bpoly.empty else np.ones_like(water)
    in_nt = contains_points(ntpoly, X, Y)
    nt_mask = water & in_bp_seg1 & in_nt
    itide[nt_mask] = 2

    dv[nt_mask, :nfmax] = value2
    dv[nt_mask, nfmax] = value1
    mfill[nt_mask] = 6
    return int(nt_mask.sum())


def set_initial_mfill(grid: MarineGrid, mfill: np.ndarray, itide: np.ndarray) -> int:
    """Set mfill=1 for tidal water only; leave land/NT untouched."""
    need = (grid.isv == 1) & (itide == 1)
    mfill[need] = 1
    return int(need.sum())



def build_node_to_triangles(triangles: np.ndarray, nnodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Build Fortran-like icells: count + cell list around each node, preserving file order."""
    counts = np.zeros(nnodes, dtype=np.int32)
    for t in range(triangles.shape[0]):
        for k in range(3):
            counts[int(triangles[t, k])] += 1
    max_count = int(counts.max()) if counts.size else 0
    adj = np.full((nnodes, max_count), -1, dtype=np.int64)
    fill = np.zeros(nnodes, dtype=np.int32)
    for t in range(triangles.shape[0]):
        for k in range(3):
            n = int(triangles[t, k])
            adj[n, fill[n]] = t
            fill[n] += 1
    return counts, adj


def compute_closest2_candidates(
    grid: MarineGrid,
    adc: AdcircData,
    need_mask: np.ndarray,
    k_query: int = 512,
) -> np.ndarray:
    """Approximate Fortran closestv2, then exact-sort candidates with row cos(lat)^2 distance.

    Fortran finds the closest two usable ADCIRC nodes for each VDatum point using
    dsq = cos(lat)^2 * dx^2 + dy^2. To avoid an impossible all-node search in Python,
    we query a generous KDTree candidate list and then exact-sort those candidates.
    """
    valid = adc.iss2.astype(bool)
    valid_idx = np.where(valid)[0]
    if valid_idx.size < 2:
        raise RuntimeError("Not enough usable ADCIRC nodes for closest-node search.")

    # Use a mid-lat scaled tree only to get candidates; final choice below uses
    # the exact Fortran row-wise c1=cos(lat)^2 metric.
    lat0 = float(np.nanmean(grid.y))
    cos0 = math.cos(RAD * lat0)
    tree_xy = np.column_stack((adc.x[valid_idx] * cos0, adc.y[valid_idx]))
    tree = cKDTree(tree_xy)

    nx, ny = grid.isv.shape
    nclo = np.full((2, nx, ny), -1, dtype=np.int64)
    k = min(k_query, valid_idx.size)

    for j, y0 in enumerate(grid.y):
        ii = np.where(need_mask[:, j])[0]
        if ii.size == 0:
            continue
        q = np.column_stack((grid.x[ii] * cos0, np.full(ii.size, y0)))
        _dist0, loc = tree.query(q, k=k)
        loc = np.atleast_2d(loc)
        if loc.shape[0] != ii.size:
            loc = loc.T
        cand = valid_idx[loc]
        c1 = math.cos(RAD * float(y0)) ** 2
        dx = adc.x[cand] - grid.x[ii][:, None]
        dy = adc.y[cand] - float(y0)
        dsq = c1 * dx * dx + dy * dy
        order = np.argsort(dsq, axis=1)
        best = np.take_along_axis(cand, order[:, :2], axis=1)
        nclo[0, ii, j] = best[:, 0]
        nclo[1, ii, j] = best[:, 1]
        if (j + 1) % 250 == 0:
            print(f"    closest2 row {j+1}/{ny}", flush=True)
    return nclo


if HAVE_NUMBA:
    @njit(cache=False)
    def _wrap_pos_angle(z):
        pi = math.atan2(0.0, -1.0)
        if z < 0.0:
            return z + 2.0 * pi
        return z

    @njit(cache=False)
    def _inside_tri3_fortran(x0, y0, x1, y1, x2, y2, x3, y3):
        # Exact translation of Fortran inside_tri3: point on side returns outside.
        a180 = math.atan2(0.0, -1.0)

        a12 = math.atan2(y2 - y1, x2 - x1)
        a13 = math.atan2(y3 - y1, x3 - x1)
        d1 = _wrap_pos_angle(a12 - a13)
        d2 = _wrap_pos_angle(a13 - a12)
        if d1 > a180:
            d = d2
            a0 = a12
        else:
            d = d1
            a0 = a13
        a1p = math.atan2(y0 - y1, x0 - x1)
        ap = _wrap_pos_angle(a1p - a0)
        if ap > d:
            return 0

        a23 = math.atan2(y3 - y2, x3 - x2)
        a21 = math.atan2(y1 - y2, x1 - x2)
        d1 = _wrap_pos_angle(a23 - a21)
        d2 = _wrap_pos_angle(a21 - a23)
        if d1 > a180:
            d = d2
            a0 = a23
        else:
            d = d1
            a0 = a21
        a2p = math.atan2(y0 - y2, x0 - x2)
        ap = _wrap_pos_angle(a2p - a0)
        if ap > d:
            return 0

        a31 = math.atan2(y1 - y3, x1 - x3)
        a32 = math.atan2(y2 - y3, x2 - x3)
        d1 = _wrap_pos_angle(a31 - a32)
        d2 = _wrap_pos_angle(a32 - a31)
        if d1 > a180:
            d = d2
            a0 = a31
        else:
            d = d1
            a0 = a32
        a3p = math.atan2(y0 - y3, x0 - x3)
        ap = _wrap_pos_angle(a3p - a0)
        if ap > d:
            return 0
        return 1

    @njit(cache=False)
    def _abc_eval_fortran(x1, y1, z1, x2, y2, z2, x3, y3, z3, xp, yp):
        # Translation of Fortran abc(): find plane z=a+b*x+c*y and evaluate.
        x = np.empty(3, dtype=np.float64)
        y = np.empty(3, dtype=np.float64)
        z = np.empty(3, dtype=np.float64)
        x[0] = x1; x[1] = x2; x[2] = x3
        y[0] = y1; y[1] = y2; y[2] = y3
        z[0] = z1; z[1] = z2; z[2] = z3
        ix = 0
        iy = 0
        for i in range(3):
            ip = i + 1
            if ip == 3:
                ip = 0
            im = i - 1
            if im == -1:
                im = 2
            if x[i] >= x[im] and x[i] <= x[ip]:
                ix = i
            if x[i] <= x[im] and x[i] >= x[ip]:
                ix = i
            if y[i] >= y[im] and y[i] <= y[ip]:
                iy = i
            if y[i] <= y[im] and y[i] >= y[ip]:
                iy = i

        im = ix - 1
        if im == -1:
            im = 2
        ip = ix + 1
        if ip == 3:
            ip = 0
        yy = y[im]
        zz1 = z[im]
        if x[ip] != x[im]:
            yy = y[im] + (x[ix] - x[im]) / (x[ip] - x[im]) * (y[ip] - y[im])
            zz1 = z[im] + (x[ix] - x[im]) / (x[ip] - x[im]) * (z[ip] - z[im])
        dzdy = (z[ix] - zz1) / (y[ix] - yy)

        im = iy - 1
        if im == -1:
            im = 2
        ip = iy + 1
        if ip == 3:
            ip = 0
        xx = x[im]
        zz2 = z[im]
        if y[ip] != y[im]:
            xx = x[im] + (y[iy] - y[im]) / (y[ip] - y[im]) * (x[ip] - x[im])
            zz2 = z[im] + (y[iy] - y[im]) / (y[ip] - y[im]) * (z[ip] - z[im])
        dzdx = (z[iy] - zz2) / (x[iy] - xx)

        a = z[0] - dzdx * x[0] - dzdy * y[0]
        return a + dzdx * xp + dzdy * yp

    @njit(cache=False)
    def _fill_interp_fortran_numba(xg, yg, node_x, node_y, triangles, ds, iss2, counts, adj, nclo, dv, mfill, nfmax):
        nx = mfill.shape[0]
        ny = mfill.shape[1]
        nfilled = 0
        for j in range(ny):
            yp = yg[j]
            for i in range(nx):
                if mfill[i, j] != 1:
                    continue
                xp = xg[i]
                did_fill = False
                for ic in range(2):
                    nc = nclo[ic, i, j]
                    if nc < 0:
                        continue
                    nall = counts[nc]
                    for nn in range(nall):
                        cell = adj[nc, nn]
                        if cell < 0:
                            continue
                        n1 = triangles[cell, 0]
                        n2 = triangles[cell, 1]
                        n3 = triangles[cell, 2]
                        zmax = ds[n1, 0]
                        if ds[n2, 0] > zmax:
                            zmax = ds[n2, 0]
                        if ds[n3, 0] > zmax:
                            zmax = ds[n3, 0]
                        if zmax > 8.999:
                            mfill[i, j] = 7
                            continue
                        if iss2[n1] + iss2[n2] + iss2[n3] < 3:
                            continue
                        inn = _inside_tri3_fortran(
                            xp, yp,
                            node_x[n1], node_y[n1],
                            node_x[n2], node_y[n2],
                            node_x[n3], node_y[n3],
                        )
                        if inn == 0:
                            continue
                        for nf in range(nfmax):
                            dv[i, j, nf] = _abc_eval_fortran(
                                node_x[n1], node_y[n1], ds[n1, nf],
                                node_x[n2], node_y[n2], ds[n2, nf],
                                node_x[n3], node_y[n3], ds[n3, nf],
                                xp, yp,
                            )
                        mfill[i, j] = 2
                        nfilled += 1
                        did_fill = True
                        break
                    if did_fill:
                        break
        return nfilled


def fill_by_triangle_interpolation(
    grid: MarineGrid,
    adc: AdcircData,
    dv: np.ndarray,
    mfill: np.ndarray,
    nfmax: int,
) -> int:
    """Fortran-faithful replacement for fill_by_interp2.

    The earlier Matplotlib TriFinder version filled the same number of cells but
    picked different triangles/values for almost every mfill=2 cell. This version
    follows vpop29nc4.f more closely: find the two closest usable ADCIRC nodes,
    then test only triangles attached to those nodes, in file order, using the
    Fortran inside_tri3 and abc plane interpolation logic.
    """
    if not HAVE_NUMBA:
        raise RuntimeError("This Fortran-style interpolation requires numba in the TAD environment.")
    need_mask = mfill == 1
    if not np.any(need_mask):
        return 0
    print("  building node-to-triangle adjacency", flush=True)
    counts, adj = build_node_to_triangles(adc.triangles, len(adc.x))
    print(f"  adjacency max triangles per node={int(counts.max())}", flush=True)
    print("  finding closest two usable ADCIRC nodes for each needed grid cell", flush=True)
    k_query = int(os.environ.get("VPOP_K_QUERY", "512"))
    print(f"  closest-node KDTree candidate count VPOP_K_QUERY={k_query}", flush=True)
    nclo = compute_closest2_candidates(grid, adc, need_mask, k_query=k_query)
    print("  *** RUNNING NUMBA FORTRAN-STYLE TRIANGLE INTERPOLATION ***", flush=True)
    return int(_fill_interp_fortran_numba(
        grid.x.astype(np.float64),
        grid.y.astype(np.float64),
        adc.x.astype(np.float64),
        adc.y.astype(np.float64),
        adc.triangles.astype(np.int64),
        adc.ds.astype(np.float64),
        adc.iss2.astype(np.int64),
        counts.astype(np.int64),
        adj.astype(np.int64),
        nclo.astype(np.int64),
        dv,
        mfill,
        int(nfmax),
    ))

def fill_small_scale(
    grid: MarineGrid,
    adc: AdcircData,
    dv: np.ndarray,
    mfill: np.ndarray,
    nfmax: int,
    max_nodes: int | None = None,
) -> int:
    need_mask = mfill == 1
    if not np.any(need_mask):
        return 0

    valid_node = adc.iss2.astype(bool) & (adc.ds[:, 0] <= 9.9)
    if not np.any(valid_node):
        return 0

    small_eps = float(os.environ.get("VPOP_SMALL_EPS", "1e-10"))

    if max_nodes is None:
        max_nodes_env = int(os.environ.get("VPOP_SMALL_MAX_NODES", "0"))
        max_nodes = max_nodes_env

    print(
        f"  small-scale settings: VPOP_SMALL_EPS={small_eps:g}, "
        f"VPOP_SMALL_MAX_NODES={max_nodes}",
        flush=True,
    )

    node_xy = np.column_stack([adc.x[valid_node], adc.y[valid_node]])
    valid_indices = np.where(valid_node)[0]
    tree = cKDTree(node_xy)

    cx = 0.5 * grid.dx
    cy = 0.5 * grid.dy

    # Use a circumscribed radius in degree space. Add eps so boundary nodes
    # that Fortran accepts are not missed by query_ball_point.
    radius = math.sqrt((cx + small_eps) * (cx + small_eps) + (cy + small_eps) * (cy + small_eps)) * 1.000001

    ii, jj = np.nonzero(need_mask)
    filled = 0

    for i, j in zip(ii, jj):
        x0 = grid.x[i]
        y0 = grid.y[j]

        cand_local = tree.query_ball_point([x0, y0], r=radius)
        if not cand_local:
            continue

        cand = valid_indices[np.asarray(cand_local, dtype=np.int64)]

        # Fortran-style node-in-cell rectangle check, with configurable tolerance.
        rect = (
            (adc.x[cand] >= x0 - cx - small_eps)
            & (adc.x[cand] <= x0 + cx + small_eps)
            & (adc.y[cand] >= y0 - cy - small_eps)
            & (adc.y[cand] <= y0 + cy + small_eps)
        )
        cand = cand[rect]
        if cand.size == 0:
            continue

        # Reproducible order. KDTree returns arbitrary candidate order.
        cand.sort()

        if max_nodes is not None and max_nodes > 0 and cand.size > max_nodes:
            cand = cand[:max_nodes]

        # Fortran local lon/lat metric: sqrt((cos(lat)*dlon)^2 + dlat^2)
        c1 = math.cos(RAD * y0) ** 2
        dist = np.sqrt(c1 * (adc.x[cand] - x0) ** 2 + (adc.y[cand] - y0) ** 2)

        zero = dist < 1.0e-14
        if np.any(zero):
            vals = adc.ds[cand[zero], :nfmax].mean(axis=0)
        else:
            w = 1.0 / dist
            vals = (adc.ds[cand, :nfmax] * w[:, None]).sum(axis=0) / w.sum()

        dv[i, j, :nfmax] = vals
        mfill[i, j] = 3
        filled += 1

    return filled


def _shift(a: np.ndarray, di: int, dj: int, fill_value=0):
    out = np.full_like(a, fill_value)
    src_i0 = max(0, -di)
    src_i1 = min(a.shape[0], a.shape[0] - di)
    dst_i0 = max(0, di)
    dst_i1 = min(a.shape[0], a.shape[0] + di)
    src_j0 = max(0, -dj)
    src_j1 = min(a.shape[1], a.shape[1] - dj)
    dst_j0 = max(0, dj)
    dst_j1 = min(a.shape[1], a.shape[1] + dj)
    out[dst_i0:dst_i1, dst_j0:dst_j1] = a[src_i0:src_i1, src_j0:src_j1]
    return out



if HAVE_NUMBA:
    @njit(cache=True)
    def _fill_weighted_fortran_sync_core(dv, mfill, itide, nfmax, wt, idx, max_iter, def_m8):
        """Faithful vpop29 fill_weighted logic.

        Matches the Fortran structure more closely than the earlier fast Step 7:
          * iade = 5,4,3,2,1
          * compute temporary tv and iupdate during one sweep
          * DO NOT let cells filled earlier in the same sweep influence later cells
          * update dv/mfill only after the sweep finishes

        Arrays use Python 0-based indexing. Interior loop 1..nx-2 / 1..ny-2
        corresponds to Fortran 2..imax2-1 / 2..jmax2-1.
        """
        fill_code = 4 if idx == 1 else 5
        nx, ny = mfill.shape
        nfp = nfmax + 1
        total = 0
        it = 0

        di_arr = np.array((-1, 1, 0, 0, -1, -1, 1, 1), dtype=np.int64)
        dj_arr = np.array((0, 0, -1, 1, -1, 1, -1, 1), dtype=np.int64)

        iupdate = np.zeros((nx, ny), dtype=np.int8)
        # Fortran tv is real*8 tv(imax2,jmax2,nfp). Keep float64.
        tv = np.empty((nx, ny, nfp), dtype=np.float64)

        for iade in range(5, 0, -1):
            while True:
                it += 1
                if it > max_iter:
                    return total, it, -999999

                # reset iupdate like Fortran: iupdate(1:imax2,1:jmax2)=0
                for i in range(nx):
                    for j in range(ny):
                        iupdate[i, j] = 0

                numit = 0
                remaining_before = 0

                # Fortran uses multisweep alternating scan directions. Because it
                # does not update mfill/dv until after the sweep, scan order should
                # not affect numerical results. Use forward order for simplicity.
                for j in range(1, ny - 1):
                    for i in range(1, nx - 1):
                        if mfill[i, j] != 1:
                            continue
                        remaining_before += 1

                        wsum1 = 0.0
                        isum3 = 0

                        for k in range(8):
                            io = i + di_arr[k]
                            jo = j + dj_arr[k]
                            mo = mfill[io, jo]
                            m1 = 0
                            if mo > 1 and mo <= 6:
                                m1 = 1
                            m2 = 0
                            if itide[io, jo] == 2:
                                m2 = 1
                            m3 = m1 * (1 - m2)
                            isum3 += m3
                            if m3 == 1:
                                wsum1 += wt[k]

                        if isum3 < iade or wsum1 <= 0.0:
                            continue

                        for nf in range(nfmax):
                            zsum1 = 0.0
                            for k in range(8):
                                io = i + di_arr[k]
                                jo = j + dj_arr[k]
                                mo = mfill[io, jo]
                                m1 = 0
                                if mo > 1 and mo <= 6:
                                    m1 = 1
                                m2 = 0
                                if itide[io, jo] == 2:
                                    m2 = 1
                                m3 = m1 * (1 - m2)
                                if m3 == 1:
                                    zsum1 += wt[k] * dv[io, jo, nf]
                            tv[i, j, nf] = zsum1 / wsum1

                        tv[i, j, nfmax] = def_m8
                        iupdate[i, j] = 1
                        numit += 1

                total += numit

                # Fortran updates dv and mfill only after whole sweep.
                if numit > 0:
                    for i in range(1, nx - 1):
                        for j in range(1, ny - 1):
                            if iupdate[i, j] == 1:
                                for nf in range(nfp):
                                    dv[i, j, nf] = tv[i, j, nf]
                                mfill[i, j] = fill_code

                remaining_after = remaining_before - numit
                print("    fortran-sync weighted iter", it, "iade", iade, "numit", numit, "numtot", total, "remaining", remaining_after)

                if remaining_before == 0:
                    return total, it, 0
                if numit == 0:
                    break

        rem = 0
        for i in range(nx):
            for j in range(ny):
                if mfill[i, j] == 1:
                    rem += 1
        return total, it, rem


def fill_weighted(
    grid: MarineGrid,
    dv: np.ndarray,
    mfill: np.ndarray,
    itide: np.ndarray,
    nfmax: int,
    idx: int,
    max_iter_factor: int = 6,
) -> int:
    
    import time

    if not HAVE_NUMBA:
        raise RuntimeError(
            "Numba is required for Fortran-sync Step 7. Install with: "
            "conda install -c conda-forge numba"
        )

    wt = get_weights(grid).astype(np.float64)
    nx, ny = mfill.shape
    max_iter = max_iter_factor * nx * ny

    if not dv.flags.c_contiguous:
        dv = np.ascontiguousarray(dv)
    if not mfill.flags.c_contiguous:
        raise RuntimeError("mfill must be C-contiguous.")
    if not itide.flags.c_contiguous:
        raise RuntimeError("itide must be C-contiguous.")

    print("    using NUMBA Fortran-SYNC weighted fill", flush=True)
    t0 = time.time()
    total, it, rem = _fill_weighted_fortran_sync_core(
        dv, mfill, itide, int(nfmax), wt, int(idx), int(max_iter), float(DEF_M8)
    )
    elapsed = time.time() - t0

    if rem == -999999:
        raise RuntimeError("Too many weighted-fill iterations; possible disconnected void.")

    print(
        f"    fortran-sync weighted fill complete: total={total}, iterations={it}, "
        f"remaining={rem}, elapsed={elapsed:.1f}s",
        flush=True,
    )
    return int(total)


def fill_layers(grid: MarineGrid, dv: np.ndarray, mfill: np.ndarray, itide: np.ndarray, nfmax: int) -> int:
    total = 0
    isvmax = int(grid.isv.max())
    if isvmax <= 1:
        return 0
    for lyr in range(2, isvmax + 1):
        layer_mask = (grid.isv == lyr) & (itide == 1)
        mfill[layer_mask] = 1
        total += fill_weighted(grid, dv, mfill, itide, nfmax, idx=2)
    return total


def write_datum_nc(filename: str | Path, grid: MarineGrid, datum: np.ndarray) -> Path:
    """Write datum NetCDF.

    Internal datum is Fortran order (imax, jmax).  Write back to NetCDF disk
    order (jmax, imax), matching the Fortran output files.
    """
    path = Path(filename)
    path = path.with_suffix(".nc")
    path.parent.mkdir(parents=True, exist_ok=True)
    imax, jmax = datum.shape
    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("box", 4)
        nc.createDimension("longitude", jmax)
        nc.createDimension("latitude", imax)
        boxv = nc.createVariable("box", "f8", ("box",))
        zv = nc.createVariable("datum", "f4", ("longitude", "latitude"), zlib=True, complevel=8)
        boxv[:] = np.array([grid.y[0], grid.x[0], grid.dy, grid.dx], dtype=np.float64)
        zv[:, :] = datum.T.astype(np.float32)
    return path


def write_byte_nc(filename: str | Path, grid: MarineGrid, field: np.ndarray, varname: str = "field") -> Path:
    """Write byte/int mask NetCDF, transposing from internal Fortran order."""
    path = Path(filename).with_suffix(".nc")
    path.parent.mkdir(parents=True, exist_ok=True)
    imax, jmax = field.shape
    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("box", 4)
        nc.createDimension("longitude", jmax)
        nc.createDimension("latitude", imax)
        boxv = nc.createVariable("box", "f8", ("box",))
        zv = nc.createVariable(varname, "i1", ("longitude", "latitude"), zlib=True, complevel=8)
        boxv[:] = np.array([grid.y[0], grid.x[0], grid.dy, grid.dx], dtype=np.float64)
        zv[:, :] = field.T.astype(np.int8)
    return path


def run_vpop(config: VpopConfig) -> dict[str, object]:
    if config.icorr_ord:
        print("NOTE: vpop29nc4.f disables icorr_ord; this Python port also skips tidal-order correction.")
    if config.icorr_dat:
        print("WARNING: station correction icorr_dat is not implemented in this vectorized draft.")
    if not _is_none(config.fdatum_msl):
        raise NotImplementedError("vpop29nc4.f stops for MSL input; Python port does the same.")

    print("STEP 1: read marine grid")
    print("  v3 mode: improved small-scale fill tolerance for mfill3/mfill4 parity", flush=True)
    grid = read_marine_grid_nc(config.vgridfile)
    print(f"  marine grid internal Fortran order: imax={grid.isv.shape[0]} jmax={grid.isv.shape[1]} isvmax={grid.isv.max()}")

    print("STEP 2: read ADCIRC datum grids")
    mode = os.environ.get("VPOP_MODE", "").strip().lower()
    if not mode:
        mode = infer_vpop_mode(config)
    x, y, tri, ds = read_adcirc_all(config.fdatum_in, config.nfmax, mode=mode)
    print(f"  ADCIRC: nodes={len(x)} triangles={len(tri)} nfmax={config.nfmax} mode={mode}")

    print("STEP 3: read polygons and screen ADCIRC nodes")
    bpoly = read_poly_segments(config.fpolygon)
    gpoly = read_poly_segments(config.gpolygon)
    ntpoly = read_poly_segments(config.nontidal)

    # Harmonize mesh longitude convention to polygon convention, like vgridder
    if not bpoly.empty:
        bp_all = np.vstack(bpoly.segments)
        bp_xmax = np.nanmax(bp_all[:, 0])

        if bp_xmax > 180.0 and np.nanmin(x) < 0.0:
            print("  wrapping ADCIRC mesh longitudes to 0..360 to match BP polygon", flush=True)
            x = wrap_lon_to_360(x)
        elif bp_xmax <= 180.0 and np.nanmax(x) > 180.0:
            print("  wrapping ADCIRC mesh longitudes to -180..180 to match BP polygon", flush=True)
            x = wrap_lon_to_180(x)

    iss, iss2 = screen_adcirc_nodes(x, y, tri, ds, bpoly, gpoly)
    adc = AdcircData(x=x, y=y, triangles=tri, ds=ds, iss=iss, iss2=iss2)
    print(f"  usable nodes: iss={int(iss.sum())} iss2={int(iss2.sum())}")

    print("STEP 4: initialize fields")
    dv, mfill, itide = initialize_fields(grid, config.nfmax)
    nt_count = mark_non_tidal(grid, bpoly, ntpoly, dv, mfill, itide, config.value1, config.value2, config.nfmax)
    ntotal = set_initial_mfill(grid, mfill, itide)
    print(f"  non-tidal cells={nt_count}; tidal water cells needing fill={ntotal}")

    debug_dir = Path(config.fdatum_out[0]).parent
    debug_dir.mkdir(parents=True, exist_ok=True)

    print("STEP 5: fill by triangle interpolation")
    n_interp = fill_by_triangle_interpolation(grid, adc, dv, mfill, config.nfmax)
    print(f"  filled by interpolation: {n_interp}; remaining={int((mfill == 1).sum())}")
    write_byte_nc(debug_dir / "mfill_after_step5.nc", grid, mfill.astype(np.int8), varname="field")

    print("STEP 6: fill small-scale cells using KDTree node-in-cell search")
    n_small = fill_small_scale(grid, adc, dv, mfill, config.nfmax)
    print(f"  filled small-scale: {n_small}; remaining={int((mfill == 1).sum())}")
    write_byte_nc(debug_dir / "mfill_after_step6.nc", grid, mfill.astype(np.int8), varname="field")

    print("STEP 7: fill remaining tidal cells by vectorized weighted neighbors")
    n_weighted = fill_weighted(grid, dv, mfill, itide, config.nfmax, idx=1)
    print(f"  filled weighted: {n_weighted}; remaining={int((mfill == 1).sum())}")
    write_byte_nc(debug_dir / "mfill_after_step7.nc", grid, mfill.astype(np.int8), varname="field")

    print("STEP 8: fill layers")
    n_layers = fill_layers(grid, dv, mfill, itide, config.nfmax)
    print(f"  filled layers: {n_layers}; remaining={int((mfill == 1).sum())}")
    write_byte_nc(debug_dir / "mfill_final_debug.nc", grid, mfill.astype(np.int8), varname="field")

    # Also write first datum snapshots for quick diagnosis.
    write_datum_nc(debug_dir / "debug_first_datum_final.nc", grid, dv[:, :, 0])

    if config.ismooth > 0:
        print("WARNING: smoothing is not implemented in this vectorized draft; skipping ismooth.")

    print("STEP 9: write outputs")
    out_paths = []
    for nf, out in enumerate(config.fdatum_out):
        out_paths.append(write_datum_nc(out, grid, dv[:, :, nf]))
    # Write itide and unfilled just like vpop29nc4.f byte outputs.
    unfilled_field = (mfill == 1).astype(np.int8)
    out_paths.append(write_byte_nc(config.unfilled, grid, unfilled_field, varname="field"))
    out_paths.append(write_byte_nc(config.fitide, grid, itide.astype(np.int8), varname="field"))
    if nt_count > 0:
        out_paths.append(write_datum_nc(config.fdatum_new, grid, dv[:, :, config.nfmax]))

    if int((mfill == 1).sum()) > 0 and config.ivoid_stop == 1:
        raise RuntimeError(f"Unfilled water cells remain: {int((mfill == 1).sum())}")

    return {
        "outputs": out_paths,
        "counts": {
            "non_tidal": nt_count,
            "need_initial": ntotal,
            "interp": n_interp,
            "small_scale": n_small,
            "weighted": n_weighted,
            "layers": n_layers,
            "remaining": int((mfill == 1).sum()),
        },
        "mfill": mfill,
        "itide": itide,
    }
def override_output_paths(cfg, out_dir: str):
    """Override all output paths to write to out_dir instead of original location."""
    os.makedirs(out_dir, exist_ok=True)
    
    # Override unfilled and itide
    cfg.unfilled = os.path.join(out_dir, "unfilled.nc")
    cfg.fitide = os.path.join(out_dir, "itide.nc")
    
    # Override datum output files
    datum_names = ["mhhw", "mhw", "mlw", "mllw", "mtl", "dtl"]
    new_fdatum_out = []
    for i in range(len(cfg.fdatum_out)):
        name = datum_names[i] if i < len(datum_names) else f"datum_{i}"
        new_fdatum_out.append(os.path.join(out_dir, f"{name}.nc"))
    cfg.fdatum_out = new_fdatum_out
    
    # Override fdatum_new if it exists
    if cfg.fdatum_new and cfg.fdatum_new.lower() != "none":
        cfg.fdatum_new = os.path.join(out_dir, os.path.basename(cfg.fdatum_new))
    
    return cfg

def main() -> int:
    print("PROGRAM vpop29 - Pacific  polygons")
    
    # Define helper function inside main
    def override_output_paths(cfg, out_dir: str):
        """Override all output paths to write to out_dir instead of original location."""
        os.makedirs(out_dir, exist_ok=True)
        
        # Override unfilled and itide
        cfg.unfilled = os.path.join(out_dir, "unfilled.nc")
        cfg.fitide = os.path.join(out_dir, "itide.nc")
        
        # Override datum output files
        datum_names = ["mhhw", "mhw", "mlw", "mllw", "mtl", "dtl"]
        new_fdatum_out = []
        for i in range(len(cfg.fdatum_out)):
            name = datum_names[i] if i < len(datum_names) else f"datum_{i}"
            new_fdatum_out.append(os.path.join(out_dir, f"{name}.nc"))
        cfg.fdatum_out = new_fdatum_out
        
        # Override fdatum_new if it exists
        if cfg.fdatum_new and cfg.fdatum_new.lower() != "none":
            cfg.fdatum_new = os.path.join(out_dir, os.path.basename(cfg.fdatum_new))
        
        return cfg
    
    # ------------------------------------------------------------
    # Get paths from environment variables
    # ------------------------------------------------------------
    vpop_indir = os.environ["VPOP_INDIR"]
    marine_grid_dir = os.environ["MARINE_GRID_DIR"]
    vpop_outdir = os.environ["VPOP_OUTDIR"]

    print(f"VPOP_INDIR: {vpop_indir}")
    print(f"MARINE_GRID_DIR: {marine_grid_dir}")
    print(f"VPOP_OUTDIR: {vpop_outdir}")
    # ------------------------------------------------------------
    # Find all polygon input files
    # ------------------------------------------------------------
    in_files = sorted(glob.glob(os.path.join(vpop_indir, "vpop29_PA*_01.in")))
    
    if not in_files:
        print(f"ERROR: No vpop29_PA*_01.in files found in {vpop_indir}")
        return 1
    
    print(f"Found {len(in_files)} polygons")
    
    # ------------------------------------------------------------
    # Get task ID from SLURM
    # ------------------------------------------------------------
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    
    if task_id < 0 or task_id >= len(in_files):
        print(f"ERROR: task_id={task_id} out of range [0, {len(in_files)-1}]")
        return 1
    
    # ------------------------------------------------------------
    # Select input file for this task
    # ------------------------------------------------------------
    infile = in_files[task_id]
    folder_name = os.path.basename(infile).replace("vpop29_", "").replace(".in", "")
    
    # ------------------------------------------------------------
    # Create output directory for this polygon
    # ------------------------------------------------------------
    out_polygon_dir = os.path.join(vpop_outdir, folder_name)
    os.makedirs(out_polygon_dir, exist_ok=True)
    
    print("=" * 50)
    print(f"SLURM Array Task ID: {task_id}")
    print(f"Polygon: {task_id + 1} of {len(in_files)}")
    print(f"Folder name: {folder_name}")
    print(f"Input file: {infile}")
    print(f"Output dir: {out_polygon_dir}")
    print("=" * 50)
    
    # ------------------------------------------------------------
    # Parse input file and override output paths
    # ------------------------------------------------------------
    cfg = parse_vpop_input(infile)
    cfg = override_output_paths(cfg, out_polygon_dir)
    
    print("Output paths:")
    print(f"  unfilled: {cfg.unfilled}")
    print(f"  itide: {cfg.fitide}")
    for i, p in enumerate(cfg.fdatum_out):
        print(f"  datum {i}: {p}")
    
    # ------------------------------------------------------------
    # Run vpop
    # ------------------------------------------------------------
    try:
        result = run_vpop(cfg)
        print("\nDONE")
        print("Counts:", result["counts"])
        print("Outputs:")
        for p in result["outputs"]:
            print(" ", p)
        return 0
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
