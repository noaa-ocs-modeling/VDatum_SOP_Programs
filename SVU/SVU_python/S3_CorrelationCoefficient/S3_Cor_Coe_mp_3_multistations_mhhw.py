#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel MATLAB-aligned MHHW correlation to multiple stations (per-shard writer; safe)

MHHW logic (station & nodes):
  • Use ONLY HHW (code = 2). In packed highs: values > 50 → decode by -100, then keep ONLY those > 50 originally.
  • Discard HW (code = 1).

Per-station MATLAB grid:
  • start = ceil(first) + 1
  • end   = floor(last) - 2
  • spacing = 20 minutes (3 per hour)
  • N = 24*3*(floor(last)-ceil(first)-3)

QC:
  • >10 events, >10 breaks, max gap ≤ 3 days → else r = 0.0 per node/station.
  • Any failure ⇒ r = 0.0 (never NaN in outputs).

Typical SLURM call:
  python S3_Cor_Coe_MHHW_matlab_aligned_mp_multiStation.py \
    --hh TADhh.nc --gage gageinfo.csv \
    --outdir out_mhhw_matlabAligned_multi \
    --out-nc coef_full_mhhw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc \
    --task-id ${SLURM_ARRAY_TASK_ID} --nchunks ${SLURM_ARRAY_TASK_COUNT} \
    --workers ${SLURM_CPUS_PER_TASK} --batch 20000 --gap-max-days 3.0 --max-stations 301
"""

import os, math, time, pathlib, argparse
import numpy as np
import xarray as xr
from scipy.interpolate import LSQUnivariateSpline
from multiprocessing import get_context

for v in ["OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(v, "1")

# ---------- Helpers (MHHW: code=2 only) ----------
def _mask_fill_array(a, fv=-9999.9):
    a = a.astype("float64")
    a[np.isclose(a, fv)] = np.nan
    return a

def _decode_highs_hhw_only(t_vec, z_vec):
    """
    Keep ONLY HHW events: packed z > 50 → decode by -100, drop others.
    Return (t_days_sorted, z_sorted, n_events, max_gap_days, span_days)
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_hhw = np.isfinite(z) & (z > 50.0)
    if not np.any(is_hhw):
        return None, None, 0, np.inf, 0.0
    z_dec = z[is_hhw] - 100.0
    t_days = (t[is_hhw] / 86400.0)
    m = np.isfinite(t_days) & np.isfinite(z_dec)
    if not np.any(m):
        return None, None, 0, np.inf, 0.0
    t_days, z_use = t_days[m], z_dec[m]
    s = np.argsort(t_days)
    t_days, z_use = t_days[s], z_use[s]
    gaps = np.diff(t_days) if t_days.size >= 2 else np.array([])
    max_gap = float(np.max(gaps)) if gaps.size else 0.0
    span = float(t_days[-1] - t_days[0]) if t_days.size >= 2 else 0.0
    return t_days, z_use, t_days.size, max_gap, span

def _fit_spline(t_days, z):
    """
    LSQ cubic spline; ~1 break per 2 days over series span.
    Return (spl_or_None, nbreaks). Caller evaluates on station grids.
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

# ---------- Station prep ----------
def _prepare_stations_MHHW(v_hh, station_idx_list, rday=1.0, np_per_hour=3, gap_max_days=3.0):
    """
    Build MATLAB-style grid and station curve for each station using MHHW (code=2 only).
    """
    K = len(station_idx_list)
    station_grids, station_curves, station_breaks = [None]*K, [None]*K, [0]*K
    stations_ok = [False]*K

    for ks, st_idx in enumerate(station_idx_list):
        row_st = v_hh.isel(station=st_idx).values  # (480,)
        t_days_st, z_st, n_ev_st, max_gap_st, _ = _decode_highs_hhw_only(row_st[:240], row_st[240:])
        if n_ev_st <= 10 or (max_gap_st > gap_max_days and (t_days_st is not None) and (t_days_st.size>1)):
            continue

        start_d = math.ceil(t_days_st[0]) + rday
        end_d   = math.floor(t_days_st[-1]) - 2*rday
        span_for_count = (math.floor(t_days_st[-1]) - math.ceil(t_days_st[0]) - 3*rday)
        n_points = int(24 * np_per_hour * max(0, span_for_count))
        if not (np.isfinite(start_d) and np.isfinite(end_d) and end_d > start_d and n_points >= 2):
            continue

        grid = np.linspace(start_d, end_d, n_points, dtype=float)
        spl_st, nbreaks_st = _fit_spline(t_days_st, z_st)
        if (spl_st is None) or (nbreaks_st <= 10):
            continue

        yy_st = spl_st(grid)
        station_grids[ks] = grid
        station_curves[ks] = yy_st
        station_breaks[ks] = nbreaks_st
        stations_ok[ks] = True

    return stations_ok, station_grids, station_curves, station_breaks

# ---------- Batch worker ----------
def _compute_block_multi(rows_batch, station_grids, station_curves, gap_max_days):
    """
    For each node in batch: decode MHHW-only, fit once, evaluate on each station grid, compute r.
    """
    B = rows_batch.shape[0]
    K = len(station_grids)
    out = np.zeros((B, K), dtype=np.float32)  # default 0.0

    for i in range(B):
        row = rows_batch[i]
        t_days, z_vals, n_ev, max_gap, _ = _decode_highs_hhw_only(row[:240], row[240:])
        if n_ev <= 10 or max_gap > gap_max_days:
            continue
        spl, nbreaks = _fit_spline(t_days, z_vals)
        if (spl is None) or (nbreaks <= 10):
            continue

        for ks, (grid_k, yy_st_k) in enumerate(zip(station_grids, station_curves)):
            if grid_k is None or yy_st_k is None:
                continue
            yy_nb_k = spl(grid_k)
            m = np.isfinite(yy_st_k) & np.isfinite(yy_nb_k)
            if not m.any():
                continue
            r = np.corrcoef(yy_st_k[m], yy_nb_k[m])[0, 1]
            out[i, ks] = np.float32(r)

    return out

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hh", required=True, help="TADhh.nc (Htime_Hval)")
    ap.add_argument("--gage", required=True, help="CSV with station_id,node_id_1b,... (use all rows up to --max-stations)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--out-nc", required=True, help="Per-shard output filename (e.g., *_SLICE_${SLURM_ARRAY_TASK_ID}.nc)")
    ap.add_argument("--batch", type=int, default=20000, help="Nodes per compute batch")
    ap.add_argument("--gap-max-days", type=float, default=3.0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--task-id", type=int)
    ap.add_argument("--nchunks", type=int)
    ap.add_argument("--max-stations", type=int, default=301)
    args = ap.parse_args()

    # stations from CSV
    import pandas as pd
    g = pd.read_csv(args.gage)
    if g.shape[0] == 0:
        raise RuntimeError("gage CSV has no rows.")
    K = min(int(args.max_stations), g.shape[0])
    station_ids_1b = g.iloc[:K, 0].astype(int).to_numpy()
    node_ids_1b    = g.iloc[:K, 1].astype(int).to_numpy()
    station_idx_list = (node_ids_1b - 1).tolist()
    print(f"[stations] using K={K} stations from CSV")

    # open highs
    ds = xr.open_dataset(args.hh, engine="netcdf4", decode_cf=False)
    v = ds["Htime_Hval"]
    if v.dims == ("high_low","station"):
        v = v.transpose("station","high_low")
    elif v.dims != ("station","high_low"):
        ds.close(); raise ValueError(f"Unexpected dims for Htime_Hval: {v.dims}")
    n_nodes_total = v.sizes["station"]
    print(f"[nc] nodes in file = {n_nodes_total:,}")

    # shard range
    task_id = args.task_id if args.task_id is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    nchunks = args.nchunks if args.nchunks is not None else int(os.environ.get("SLURM_ARRAY_TASK_COUNT", "1"))
    nodes_per = math.ceil(n_nodes_total / max(1, nchunks))
    i0 = task_id * nodes_per
    i1 = min(n_nodes_total, (task_id + 1) * nodes_per)
    if i1 <= i0:
        ds.close(); print(f"[shard] empty shard {task_id}/{nchunks}"); return
    print(f"[shard] {task_id+1}/{nchunks}: node range [{i0}:{i1})")

    # station prep
    print("[prep] building station grids & curves (MHHW, MATLAB-aligned)…")
    stations_ok, station_grids, station_curves, station_breaks = _prepare_stations_MHHW(
        v, station_idx_list, rday=1.0, np_per_hour=3, gap_max_days=args.gap_max_days
    )
    print(f"[prep] stations OK: {sum(stations_ok)}/{K}")

    # compute shard
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
                rblk = np.vstack(parts).astype(np.float32, copy=False)
            else:
                rblk = _compute_block_multi(rows, station_grids, station_curves, args.gap_max_days)

            coef_shard[(idx_batch_nodes - i0).astype(np.int64), :] = rblk

            processed += idx_batch_nodes.size
            now = time.time()
            if (now - last_log) >= 2.0 or bstop == idx_all.size:
                elapsed = (now - t_start) / 60.0
                pct = 100.0 * processed / idx_all.size
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

    # save shard file
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    out_path = str(pathlib.Path(args.outdir) / args.out_nc)
    xr.Dataset(
        data_vars=dict(coef=(("node","station"), coef_shard)),
        coords=dict(
            node=("node", np.arange(i0, i1, dtype=np.int64)),
            station=("station", station_ids_1b.astype(np.int64))
        ),
        attrs=dict(datum="mhhw", ns=i0, ne=i1, note="Per-shard MHHW coefficients; MATLAB-aligned; multi-station")
    ).to_netcdf(out_path, mode="w", encoding={"coef":{"dtype":"float32"}})
    print(f"[save] wrote shard file: {out_path}")

if __name__ == "__main__":
    main()

