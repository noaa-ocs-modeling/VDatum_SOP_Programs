#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parallel MATLAB-aligned MTL / DTL correlation to multiple stations (per-shard writer; safe)

Datums:

  MTL:
    - High family: HW + HHW  (codes 1 or 2)
    - Low  family: HLW + LLW (codes -1 or -2)
    - Curve = (spline_high + spline_low) / 2 on a common grid spanning both families.

  DTL:
    - High family: HHW  only (code  2)
    - Low  family: LLW  only (code -2)
    - Curve = (spline_HHW + spline_LLW) / 2 on a common grid spanning both.

Both families use packed high/low series from:
  - Htime_Hval (TADhh.nc)
  - Ltime_Lval (TADll.nc)

Per-station MATLAB-like grid:
  - start = ceil(min(t_high, t_low)) + rday
  - end   = floor(max(t_high, t_low)) - 2*rday
  - spacing = 20 min (np_per_hour=3)
  - N = 24*np_per_hour*(floor(tmax)-ceil(tmin)-3*rday)

QC:
  - Each family (high, low) must have >10 events, max_gap <= gap_max_days, spline with >10 breaks.
  - Any failure → curve becomes NaN → correlation r = 0.0 for that node/station.

Output (per shard):
  - NetCDF with coef(node_in_shard, station=K), plus node and station coords.

Call example:

  python S3_Cor_Coe_MTL_DTL_matlab_aligned_mp_multiStation.py \
      --hh TADhh.nc --ll TADll.nc --gage gageinfo.csv \
      --datum mtl \
      --outdir out_mtl_matlabAligned_multi \
      --out-nc coef_full_mtl_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc \
      --task-id ${SLURM_ARRAY_TASK_ID} --nchunks ${SLURM_ARRAY_TASK_COUNT} \
      --workers ${SLURM_CPUS_PER_TASK} --batch 20000 --gap-max-days 3.0 --max-stations 301

Then repeat with --datum dtl and a different outdir/out-nc.
"""
import os, math, time, pathlib, argparse
import numpy as np
import xarray as xr
from scipy.interpolate import LSQUnivariateSpline
from multiprocessing import get_context

# single-thread math per worker
for v in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ.setdefault(v, "1")


# ---------------- basic helpers ----------------

def _mask_fill_array(a, fv=-9999.9):
    """
    Convert fill values to NaN and ensure float64.
    """
    a = a.astype("float64")
    a[np.isclose(a, fv)] = np.nan
    return a


def _decode_highs_all(t_vec, z_vec):
    """
    Decode highs: HW + HHW (MHW family).
      - packed HHW: z > 50 → z_dec = z - 100
      - packed HW:  all other finite z
    Use ALL finite highs after decoding.

    Returns:
      t_days_sorted, z_sorted, n_events, max_gap_days, span_days
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_hhw = np.isfinite(z) & (z > 50.0)
    z_dec = np.where(is_hhw, z - 100.0, z)
    m = np.isfinite(t) & np.isfinite(z_dec)
    if not np.any(m):
        return None, None, 0, np.inf, 0.0
    t_days = t[m] / 86400.0
    z_use = z_dec[m]
    s = np.argsort(t_days)
    t_days, z_use = t_days[s], z_use[s]
    gaps = np.diff(t_days) if t_days.size >= 2 else np.array([])
    max_gap = float(np.max(gaps)) if gaps.size else 0.0
    span = float(t_days[-1] - t_days[0]) if t_days.size >= 2 else 0.0
    return t_days, z_use, t_days.size, max_gap, span


def _decode_highs_hhw_only(t_vec, z_vec):
    """
    Decode highs: HHW only (MHHW family).
      - packed HHW: z > 50 → z_dec = z - 100
      - ignore all others
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_hhw = np.isfinite(z) & (z > 50.0)
    if not np.any(is_hhw):
        return None, None, 0, np.inf, 0.0
    t_days = t[is_hhw] / 86400.0
    z_dec = z[is_hhw] - 100.0
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


def _decode_lows_all(t_vec, z_vec):
    """
    Decode lows: HLW + LLW (MLW family).
      - packed LLW: z < -50 → z_dec = z + 100
      - packed HLW: all other finite z
    Use ALL finite lows after decoding.
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_llw = np.isfinite(z) & (z < -50.0)
    z_dec = np.where(is_llw, z + 100.0, z)
    m = np.isfinite(t) & np.isfinite(z_dec)
    if not np.any(m):
        return None, None, 0, np.inf, 0.0
    t_days = t[m] / 86400.0
    z_use = z_dec[m]
    s = np.argsort(t_days)
    t_days, z_use = t_days[s], z_use[s]
    gaps = np.diff(t_days) if t_days.size >= 2 else np.array([])
    max_gap = float(np.max(gaps)) if gaps.size else 0.0
    span = float(t_days[-1] - t_days[0]) if t_days.size >= 2 else 0.0
    return t_days, z_use, t_days.size, max_gap, span


def _decode_lows_llw_only(t_vec, z_vec):
    """
    Decode lows: LLW only (MLLW family).
      - packed LLW: z < -50 → z_dec = z + 100
      - ignore others
    """
    t = _mask_fill_array(t_vec)
    z = _mask_fill_array(z_vec)
    is_llw = np.isfinite(z) & (z < -50.0)
    if not np.any(is_llw):
        return None, None, 0, np.inf, 0.0
    t_days = t[is_llw] / 86400.0
    z_dec = z[is_llw] + 100.0
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
    LSQ cubic spline; ~1 break per 2 days over the series span.
    Returns (spl_or_None, nbreaks).
    """
    if t_days is None or z is None or t_days.size == 0:
        return None, 0
    span = float(t_days[-1] - t_days[0])
    nbreaks = max(2, int(round(span / 2.0)))  # ~1 break per 2 days
    knots = []
    if nbreaks > 2:
        br = np.linspace(t_days[0], t_days[-1], nbreaks)
        knots = br[1:-1]
    try:
        spl = LSQUnivariateSpline(t_days, z, t=knots, k=3)
        return spl, nbreaks
    except Exception:
        return None, nbreaks


# ---------------- station prep (MTL / DTL) ----------------

def _prepare_stations_MTL_DTL(v_hh, v_ll, station_idx_list, datum, rday=1.0,
                              np_per_hour=3, gap_max_days=3.0):
    """
    Build MATLAB-style grid and station curves for each station:

      datum = "mtl":
        high decoder = _decode_highs_all
        low  decoder = _decode_lows_all

      datum = "dtl":
        high decoder = _decode_highs_hhw_only
        low  decoder = _decode_lows_llw_only

    For each station:
      - Build a common time grid spanning both families.
      - Fit a spline to each family.
      - Station curve = average of families on that grid.
    """
    assert datum in ("mtl", "dtl")
    K = len(station_idx_list)
    station_grids = [None] * K
    station_curves = [None] * K
    station_breaks_hi = [0] * K
    station_breaks_lo = [0] * K
    stations_ok = [False] * K

    for ks, st_idx in enumerate(station_idx_list):
        row_h = v_hh.isel(station=st_idx).values  # (480,)
        row_l = v_ll.isel(station=st_idx).values  # (480,)

        if datum == "mtl":
            th, zh, n_ev_h, max_gap_h, _ = _decode_highs_all(row_h[:240], row_h[240:])
            tl, zl, n_ev_l, max_gap_l, _ = _decode_lows_all(row_l[:240], row_l[240:])
        else:  # dtl
            th, zh, n_ev_h, max_gap_h, _ = _decode_highs_hhw_only(row_h[:240], row_h[240:])
            tl, zl, n_ev_l, max_gap_l, _ = _decode_lows_llw_only(row_l[:240], row_l[240:])

        # require both families to have enough events and not huge gaps
        if (n_ev_h <= 10) or (n_ev_l <= 10):
            continue
        if (max_gap_h > gap_max_days and th is not None and th.size > 1):
            continue
        if (max_gap_l > gap_max_days and tl is not None and tl.size > 1):
            continue

        # Build common grid spanning both families
        t_concat = []
        if th is not None and th.size > 0:
            t_concat.append(th)
        if tl is not None and tl.size > 0:
            t_concat.append(tl)
        if not t_concat:
            continue
        t_all = np.concatenate(t_concat)
        tmin = float(np.nanmin(t_all))
        tmax = float(np.nanmax(t_all))

        start_d = math.ceil(tmin) + rday
        end_d = math.floor(tmax) - 2.0 * rday
        span_for_count = (math.floor(tmax) - math.ceil(tmin) - 3.0 * rday)
        n_points = int(24 * np_per_hour * max(0.0, span_for_count))

        if not (np.isfinite(start_d) and np.isfinite(end_d) and end_d > start_d and n_points >= 2):
            continue

        grid = np.linspace(start_d, end_d, n_points, dtype=float)

        # Fit station high family
        spl_h, nbreaks_h = _fit_spline(th, zh)
        if (spl_h is None) or (nbreaks_h <= 10):
            continue
        yy_h = spl_h(grid)

        # Fit station low family
        spl_l, nbreaks_l = _fit_spline(tl, zl)
        if (spl_l is None) or (nbreaks_l <= 10):
            continue
        yy_l = spl_l(grid)

        yy = 0.5 * (yy_h + yy_l)

        station_grids[ks] = grid
        station_curves[ks] = yy
        station_breaks_hi[ks] = nbreaks_h
        station_breaks_lo[ks] = nbreaks_l
        stations_ok[ks] = True

    return stations_ok, station_grids, station_curves, station_breaks_hi, station_breaks_lo


# ---------------- batch worker ----------------

def _decode_node_families(row_h, row_l, datum, gap_max_days):
    """
    Decode high and low families for a single node according to datum.
    Returns:
      (th, zh, n_ev_h, max_gap_h, tl, zl, n_ev_l, max_gap_l, ok)
    where ok includes the basic QC on gaps, but not spline fitting.
    """
    if datum == "mtl":
        th, zh, n_ev_h, max_gap_h, _ = _decode_highs_all(row_h[:240], row_h[240:])
        tl, zl, n_ev_l, max_gap_l, _ = _decode_lows_all(row_l[:240], row_l[240:])
    else:  # dtl
        th, zh, n_ev_h, max_gap_h, _ = _decode_highs_hhw_only(row_h[:240], row_h[240:])
        tl, zl, n_ev_l, max_gap_l, _ = _decode_lows_llw_only(row_l[:240], row_l[240:])

    if (n_ev_h <= 10) or (n_ev_l <= 10):
        return None, None, 0, 0.0, None, None, 0, 0.0, False
    if (max_gap_h > gap_max_days and th is not None and th.size > 1):
        return None, None, 0, 0.0, None, None, 0, 0.0, False
    if (max_gap_l > gap_max_days and tl is not None and tl.size > 1):
        return None, None, 0, 0.0, None, None, 0, 0.0, False

    return th, zh, n_ev_h, max_gap_h, tl, zl, n_ev_l, max_gap_l, True


def _compute_block_multi_MTL_DTL(rows_hh_batch, rows_ll_batch,
                                 station_grids, station_curves,
                                 gap_max_days, datum):
    """
    For each node in the batch:
      - decode high and low families
      - fit two splines
      - evaluate on each station's grid
      - average families and compute Pearson r
    """
    B = rows_hh_batch.shape[0]
    K = len(station_grids)
    out = np.zeros((B, K), dtype=np.float32)  # default 0.0

    for i in range(B):
        row_h = rows_hh_batch[i]
        row_l = rows_ll_batch[i]

        th, zh, n_ev_h, _, tl, zl, n_ev_l, _, ok = _decode_node_families(row_h, row_l, datum, gap_max_days)
        if not ok:
            continue

        spl_h, nbreaks_h = _fit_spline(th, zh)
        if (spl_h is None) or (nbreaks_h <= 10):
            continue
        spl_l, nbreaks_l = _fit_spline(tl, zl)
        if (spl_l is None) or (nbreaks_l <= 10):
            continue

        for ks, (grid_k, yy_st_k) in enumerate(zip(station_grids, station_curves)):
            if grid_k is None or yy_st_k is None:
                continue
            yy_h = spl_h(grid_k)
            yy_l = spl_l(grid_k)
            yy_nb_k = 0.5 * (yy_h + yy_l)
            m = np.isfinite(yy_st_k) & np.isfinite(yy_nb_k)
            if not m.any():
                continue
            r = np.corrcoef(yy_st_k[m], yy_nb_k[m])[0, 1]
            if np.isfinite(r):
                out[i, ks] = np.float32(r)

    return out


# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hh", required=True, help="TADhh.nc (Htime_Hval)")
    ap.add_argument("--ll", required=True, help="TADll.nc (Ltime_Lval)")
    ap.add_argument("--gage", required=True, help="CSV: station_id,node_id_1b,...; first max-stations used")
    ap.add_argument("--datum", required=True, choices=["mtl", "dtl"],
                    help="Which datum to compute ('mtl' or 'dtl')")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--out-nc", required=True, help="Per-shard output filename")
    ap.add_argument("--batch", type=int, default=20000, help="Nodes per compute batch")
    ap.add_argument("--gap-max-days", type=float, default=3.0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--task-id", type=int)
    ap.add_argument("--nchunks", type=int)
    ap.add_argument("--max-stations", type=int, default=301)
    args = ap.parse_args()

    datum = args.datum.lower()

    # stations from CSV
    import pandas as pd
    g = pd.read_csv(args.gage)
    if g.shape[0] == 0:
        raise RuntimeError("gage CSV has no rows.")
    K = min(int(args.max_stations), g.shape[0])
    station_ids_1b = g.iloc[:K, 0].astype(int).to_numpy()
    node_ids_1b = g.iloc[:K, 1].astype(int).to_numpy()
    station_idx_list = (node_ids_1b - 1).tolist()
    print(f"[stations] using K={K} stations from CSV")

    # open highs & lows
    ds_hh = xr.open_dataset(args.hh, engine="netcdf4", decode_cf=False)
    ds_ll = xr.open_dataset(args.ll, engine="netcdf4", decode_cf=False)

    v_hh = ds_hh["Htime_Hval"]
    v_ll = ds_ll["Ltime_Lval"]

    if v_hh.dims == ("high_low", "station"):
        v_hh = v_hh.transpose("station", "high_low")
    elif v_hh.dims != ("station", "high_low"):
        ds_hh.close(); ds_ll.close()
        raise ValueError(f"Unexpected dims for Htime_Hval: {v_hh.dims}")

    if v_ll.dims == ("high_low", "station"):
        v_ll = v_ll.transpose("station", "high_low")
    elif v_ll.dims != ("station", "high_low"):
        ds_hh.close(); ds_ll.close()
        raise ValueError(f"Unexpected dims for Ltime_Lval: {v_ll.dims}")

    n_nodes_total = v_hh.sizes["station"]
    if n_nodes_total != v_ll.sizes["station"]:
        ds_hh.close(); ds_ll.close()
        raise RuntimeError("Htime_Hval and Ltime_Lval have different station sizes.")
    print(f"[nc] nodes in file = {n_nodes_total:,}")

    # shard range
    task_id = args.task_id if args.task_id is not None else int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    nchunks = args.nchunks if args.nchunks is not None else int(os.environ.get("SLURM_ARRAY_TASK_COUNT", "1"))
    nodes_per = math.ceil(n_nodes_total / max(1, nchunks))
    i0 = task_id * nodes_per
    i1 = min(n_nodes_total, (task_id + 1) * nodes_per)
    if i1 <= i0:
        ds_hh.close(); ds_ll.close()
        print(f"[shard] empty shard {task_id}/{nchunks}")
        return
    print(f"[shard] {task_id+1}/{nchunks}: node range [{i0}:{i1})")

    # station prep
    print(f"[prep] building station grids & curves ({datum.upper()}, MATLAB-aligned)…")
    stations_ok, station_grids, station_curves, br_hi, br_lo = _prepare_stations_MTL_DTL(
        v_hh, v_ll, station_idx_list,
        datum=datum, rday=1.0, np_per_hour=3, gap_max_days=args.gap_max_days
    )
    print(f"[prep] stations OK: {sum(stations_ok)}/{K}")

    coef_shard = np.zeros((i1 - i0, K), dtype=np.float32)  # default 0.0
    BATCH = int(args.batch)
    WORKERS = max(1, int(args.workers))
    t_start = time.time(); processed = 0; last_log = time.time()

    def read_rows(idx_list):
        rows_hh = v_hh.isel(station=xr.DataArray(idx_list, dims="sidx")).values
        rows_ll = v_ll.isel(station=xr.DataArray(idx_list, dims="sidx")).values
        return rows_hh, rows_ll

    pool = None
    if WORKERS > 1:
        ctx = get_context("fork")
        pool = ctx.Pool(processes=WORKERS)

    try:
        idx_all = np.arange(i0, i1, dtype=np.int64)
        for bstart in range(0, idx_all.size, BATCH):
            bstop = min(bstart + BATCH, idx_all.size)
            idx_batch_nodes = idx_all[bstart:bstop]
            rows_hh_batch, rows_ll_batch = read_rows(idx_batch_nodes)

            if WORKERS > 1:
                chunks_hh = np.array_split(rows_hh_batch, WORKERS, axis=0)
                chunks_ll = np.array_split(rows_ll_batch, WORKERS, axis=0)
                work = [
                    (chh, cll, station_grids, station_curves, args.gap_max_days, datum)
                    for chh, cll in zip(chunks_hh, chunks_ll)
                ]
                parts = pool.starmap(_compute_block_multi_MTL_DTL, work)
                rblk = np.vstack(parts).astype(np.float32, copy=False)
            else:
                rblk = _compute_block_multi_MTL_DTL(
                    rows_hh_batch, rows_ll_batch,
                    station_grids, station_curves,
                    args.gap_max_days, datum
                )

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
        ds_hh.close(); ds_ll.close()

    # save shard file
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    out_path = str(pathlib.Path(args.outdir) / args.out_nc)
    xr.Dataset(
        data_vars=dict(coef=(("node", "station"), coef_shard)),
        coords=dict(
            node=("node", np.arange(i0, i1, dtype=np.int64)),
            station=("station", station_ids_1b.astype(np.int64)),
        ),
        attrs=dict(datum=datum, ns=i0, ne=i1,
                   note=f"Per-shard {datum.upper()} coefficients; MATLAB-aligned; multi-station")
    ).to_netcdf(out_path, mode="w", encoding={"coef": {"dtype": "float32"}})
    print(f"[save] wrote shard file: {out_path}")


if __name__ == "__main__":
    main()


