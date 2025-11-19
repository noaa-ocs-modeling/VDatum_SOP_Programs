#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel MATLAB-aligned MLW correlation to *multiple stations* (per-shard writer; safe)

Exact math per station (MLW):
  • Decode LLW (+100), then use ALL finite lows (MLW + MLLW) by code (−1 or −2).
  • Station-defined grid for *that station*:
       start = ceil(first) + 1
       end   = floor(last) - 2
       spacing = 20 minutes (3 per hour)
       N = 24*3*(floor(last)-ceil(first)-3)
  • Node QC: >10 events, >10 breaks, max gap ≤ 3 days → else r = 0.0 for that node/station.
  • Any failure ⇒ r = 0.0 (never NaN in outputs).

Parallelization:
  • Shard nodes across SLURM array (task-id / nchunks).
  • Optional local multiprocessing over node batches.
  • Each shard writes ONE NetCDF: coef(node_in_shard, station=K), coords include absolute node ids and station IDs.

Usage (typical SLURM srun):
  python S3_Cor_Coe_MLW_matlab_aligned_mp_multiStation.py \
    --ll TADll.nc --gage gageinfo.csv \
    --outdir out_mlw_matlabAligned \
    --out-nc coef_full_mlw_matlabAligned_SLICE_${SLURM_ARRAY_TASK_ID}.nc \
    --task-id ${SLURM_ARRAY_TASK_ID} --nchunks ${SLURM_ARRAY_TASK_COUNT} \
    --workers ${SLURM_CPUS_PER_TASK} --batch 20000 --gap-max-days 3.0 --max-stations 301
"""

import os, math, time, pathlib, argparse
import numpy as np
import xarray as xr
from scipy.interpolate import LSQUnivariateSpline
from multiprocessing import get_context

# single-thread math per worker
for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(v, "1")

# ---------------- Helpers (unchanged MLW math) ----------------
def _mask_fill_array(a, fv=-9999.9):
    a = a.astype("float64")
    a[np.isclose(a, fv)] = np.nan
    return a

def _decode_lows_use_all(t_vec, z_vec):
    """
    240 times (sec) + 240 lows (packed). Decode LLW by +100 (z<-50),
    then use ALL finite lows (MLW+MLLW), sort by time (days), and report stats.
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_llw = np.isfinite(z) & (z < -50.0)
    z_dec  = np.where(is_llw, z + 100.0, z)
    m = np.isfinite(t) & np.isfinite(z_dec)
    if not np.any(m):
        return None, None, 0, np.inf, 0.0
    t_days = t[m] / 86400.0
    z_use  = z_dec[m]
    s = np.argsort(t_days)
    t_days, z_use = t_days[s], z_use[s]
    gaps = np.diff(t_days) if t_days.size >= 2 else np.array([])
    max_gap = float(np.max(gaps)) if gaps.size else 0.0
    span = float(t_days[-1] - t_days[0]) if t_days.size >= 2 else 0.0
    return t_days, z_use, t_days.size, max_gap, span

def _fit_spline(t_days, z):
    """
    Fit LSQ cubic spline; ~1 break per 2 days over node's own span.
    Return (spl_or_None, nbreaks). Caller will evaluate on station-specific grids.
    """
    if t_days is None or z is None or t_days.size == 0:
        return None, 0
    span = float(t_days[-1] - t_days[0])
    nbreaks = max(2, int(round(span / 2.0)))
    knots = []
    if nbreaks > 2:
        br = np.linspace(t_days[0], t_days[-1], nbreaks)
        knots = br[1:-1]
    try:
        spl = LSQUnivariateSpline(t_days, z, t=knots, k=3)
        return spl, nbreaks
    except Exception:
        return None, nbreaks

# ---------------- Station prep ----------------
def _prepare_stations_MLW(v_ll, station_idx_list, rday=1.0, np_per_hour=3, gap_max_days=3.0):
    """
    For each station index, build its MATLAB-style grid and station spline curve.
    Returns:
      stations_ok   : list of booleans (station passes QC)
      station_ids   : np.array shape (K,)
      station_grids : list of 1D arrays (grid per station; or None if failed)
      station_curves: list of 1D arrays (yy_st per station on its own grid; or None)
      station_breaks: list of integers (nbreaks per station; 0 if failed)
    """
    K = len(station_idx_list)
    stations_ok = [False]*K
    station_grids, station_curves, station_breaks = [None]*K, [None]*K, [0]*K

    for ks, st_idx in enumerate(station_idx_list):
        row_st = v_ll.isel(station=st_idx).values  # (480,)
        t_days_st, z_st, n_ev_st, max_gap_st, span_st = _decode_lows_use_all(row_st[:240], row_st[240:])
        if n_ev_st <= 10 or (max_gap_st > gap_max_days and (t_days_st is not None) and (t_days_st.size>1)):
            # station fails QC => mark as bad, keep None grid/curve
            stations_ok[ks] = False
            continue

        start_d = math.ceil(t_days_st[0]) + rday
        end_d   = math.floor(t_days_st[-1]) - 2*rday
        span_for_count = (math.floor(t_days_st[-1]) - math.ceil(t_days_st[0]) - 3*rday)
        n_points = int(24 * np_per_hour * max(0, span_for_count))
        if not (np.isfinite(start_d) and np.isfinite(end_d) and end_d > start_d and n_points >= 2):
            stations_ok[ks] = False
            continue

        grid = np.linspace(start_d, end_d, n_points, dtype=float)
        spl_st, nbreaks_st = _fit_spline(t_days_st, z_st)
        if (spl_st is None) or (nbreaks_st <= 10):
            stations_ok[ks] = False
            continue

        yy_st = spl_st(grid)  # evaluate on its own grid
        station_grids[ks] = grid
        station_curves[ks] = yy_st
        station_breaks[ks] = nbreaks_st
        stations_ok[ks] = True

    return stations_ok, station_grids, station_curves, station_breaks

# ---------------- Batch worker ----------------
def _compute_block_multi(rows_batch, station_grids, station_curves, gap_max_days):
    """
    Compute r for a batch of node rows against ALL prepared stations.
    Returns r_batch: shape (B, K)
    """
    B = rows_batch.shape[0]
    K = len(station_grids)
    out = np.zeros((B, K), dtype=np.float32)  # default 0.0

    # For each node (once): decode, QC, fit ONE spline, then evaluate on each station grid
    for i in range(B):
        row = rows_batch[i]
        t_days, z_vals, n_ev, max_gap, span = _decode_lows_use_all(row[:240], row[240:])
        if n_ev <= 10 or max_gap > gap_max_days:
            continue
        spl, nbreaks = _fit_spline(t_days, z_vals)
        if (spl is None) or (nbreaks <= 10):
            continue

        for ks, (grid_k, yy_st_k) in enumerate(zip(station_grids, station_curves)):
            if grid_k is None or yy_st_k is None:
                # station failed prep → leave 0.0
                continue
            yy_nb_k = spl(grid_k)  # evaluate node spline on station's grid
            m = np.isfinite(yy_st_k) & np.isfinite(yy_nb_k)
            if not m.any():
                continue
            r = np.corrcoef(yy_st_k[m], yy_nb_k[m])[0, 1]
            out[i, ks] = np.float32(r)

    return out

# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ll", required=True, help="TADll.nc (Ltime_Lval)")
    ap.add_argument("--gage", required=True, help="CSV with station_id,node_id,... (use all rows up to --max-stations)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--out-nc", required=True, help="Per-shard output filename (e.g., *_SLICE_${SLURM_ARRAY_TASK_ID}.nc)")
    ap.add_argument("--batch", type=int, default=20000, help="Nodes per compute batch")
    ap.add_argument("--gap-max-days", type=float, default=3.0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--task-id", type=int)
    ap.add_argument("--nchunks", type=int)
    ap.add_argument("--max-stations", type=int, default=301)
    args = ap.parse_args()

    # ---- Read stations from CSV (station_id, node_id_1b in first 2 columns) ----
    import pandas as pd
    g = pd.read_csv(args.gage)
    if g.shape[0] == 0:
        raise RuntimeError("gage CSV has no rows.")
    K = min(int(args.max_stations), g.shape[0])
    station_ids_1b = g.iloc[:K, 0].astype(int).to_numpy()
    node_ids_1b    = g.iloc[:K, 1].astype(int).to_numpy()
    station_idx_list = (node_ids_1b - 1).tolist()
    print(f"[stations] using K={K} stations from CSV")

    # ---- Open lows once ----
    ds = xr.open_dataset(args.ll, engine="netcdf4", decode_cf=False)
    v = ds["Ltime_Lval"]
    if v.dims == ("high_low","station"):
        v = v.transpose("station","high_low")
    elif v.dims != ("station","high_low"):
        ds.close(); raise ValueError(f"Unexpected dims for Ltime_Lval: {v.dims}")
    n_nodes_total = v.sizes["station"]
    print(f"[nc] nodes in file = {n_nodes_total:,}")

    # ---- Shard range ----
    task_id = args.task_id if args.task_id is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    nchunks = args.nchunks if args.nchunks is not None else int(os.environ.get("SLURM_ARRAY_TASK_COUNT", "1"))
    nodes_per = math.ceil(n_nodes_total / max(1, nchunks))
    i0 = task_id * nodes_per
    i1 = min(n_nodes_total, (task_id + 1) * nodes_per)
    if i1 <= i0:
        ds.close(); print(f"[shard] empty shard {task_id}/{nchunks}"); return
    print(f"[shard] {task_id+1}/{nchunks}: node range [{i0}:{i1})")

    # ---- Prepare each station's grid and station spline curve ----
    print("[prep] building station grids & curves (MLW, MATLAB-aligned)…")
    stations_ok, station_grids, station_curves, station_breaks = _prepare_stations_MLW(
        v, station_idx_list, rday=1.0, np_per_hour=3, gap_max_days=args.gap_max_days
    )
    ok_count = sum(stations_ok)
    print(f"[prep] stations OK: {ok_count}/{K}")

    # ---- Compute this shard → coef_shard (i1-i0, K) ----
    coef_shard = np.zeros((i1 - i0, K), dtype=np.float32)  # default 0.0
    BATCH = int(args.batch)
    WORKERS = max(1, int(args.workers))
    t_start = time.time(); processed = 0; last_log = time.time()

    def read_rows(idx_list):
        return v.isel(station=xr.DataArray(idx_list, dims="sidx")).values  # (batch, 480)

    pool = None
    if WORKERS > 1:
        ctx = get_context("fork")
        pool = ctx.Pool(processes=WORKERS)

    try:
        idx_all = np.arange(i0, i1, dtype=np.int64)
        for bstart in range(0, idx_all.size, BATCH):
            bstop = min(bstart + BATCH, idx_all.size)
            idx_batch_nodes = idx_all[bstart:bstop]
            rows = read_rows(idx_batch_nodes)

            if WORKERS > 1:
                chunks = np.array_split(rows, WORKERS, axis=0)
                parts = pool.starmap(_compute_block_multi, [(ch, station_grids, station_curves, args.gap_max_days) for ch in chunks])
                rblk = np.vstack(parts).astype(np.float32, copy=False)  # (batch, K)
            else:
                rblk = _compute_block_multi(rows, station_grids, station_curves, args.gap_max_days)

            rel = (idx_batch_nodes - i0).astype(np.int64)
            coef_shard[rel, :] = rblk

            processed += idx_batch_nodes.size
            now = time.time()
            if (now - last_log) >= 2.0 or bstop == idx_all.size:
                elapsed = (now - t_start) / 60.0
                pct = 100.0 * processed / idx_all.size
                # quick stats over this r-block
                finite = np.isfinite(rblk)
                rmin = float(np.min(rblk[finite])) if finite.any() else float("nan")
                rmax = float(np.max(rblk[finite])) if finite.any() else float("nan")
                rmean = float(np.mean(rblk[finite])) if finite.any() else float("nan")
                n_zero = int(np.count_nonzero(rblk == 0.0))
                print(f"[batch] {bstart}:{bstop}  progress={pct:5.1f}%  elapsed={elapsed:7.2f} min  "
                      f"zeros={n_zero:7d}  r(min/mean/max)=({rmin:.3f},{rmean:.3f},{rmax:.3f})")
                last_log = now
    finally:
        if pool is not None:
            pool.close(); pool.join()
        ds.close()

    # ---- Save ONE file for this shard (safe) ----
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    out_path = str(pathlib.Path(args.outdir) / args.out_nc)
    xr.Dataset(
        data_vars=dict(coef=(("node","station"), coef_shard)),
        coords=dict(
            node=("node", np.arange(i0, i1, dtype=np.int64)),          # absolute 0-based node ids
            station=("station", station_ids_1b.astype(np.int64))       # station IDs from CSV
        ),
        attrs=dict(datum="mlw", ns=i0, ne=i1, note="Per-shard MLW coefficients; MATLAB-aligned; multi-station")
    ).to_netcdf(out_path, mode="w", encoding={"coef":{"dtype":"float32"}})
    print(f"[save] wrote shard file: {out_path}")
if __name__ == "__main__":
    main()





