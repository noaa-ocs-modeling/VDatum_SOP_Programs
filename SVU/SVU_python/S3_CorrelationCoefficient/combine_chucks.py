#!/usr/bin/env python3
"""
Combine per-shard MlW NetCDFs into a single (node, station) file.

Assumes shard files look like:
    out_mlw_matlabAligned_multi/coef_full_mlw_matlabAligned_multi_SLICE_*.nc

Output:
    coef_full_mlw_matlabAligned_multi.nc
"""

import os
import glob
import numpy as np
import xarray as xr


def main():
    outdir = "out_mhhw_matlabAligned_multi"
    pattern = "coef_full_mhhw_matlabAligned_multi_SLICE_*.nc"
    parts = sorted(glob.glob(os.path.join(outdir, pattern)))

    if not parts:
        raise RuntimeError(f"No shard files found matching {os.path.join(outdir, pattern)}")

    # Get station ids from the first shard and assert consistency across shards
    with xr.open_dataset(parts[0]) as d0:
        station_ids = d0["station"].values.astype(np.int64)
        K = station_ids.size
    print(f"[combine] Found {len(parts)} shard files; K={K} stations")

    # Discover total number of nodes (assuming node coord is 0-based index)
    max_node = -1
    for p in parts:
        with xr.open_dataset(p) as ds:
            nmax = int(ds["node"].values.max())
            if nmax > max_node:
                max_node = nmax
    N = max_node + 1
    print(f"[combine] Total nodes inferred: N={N}")

    # Assemble full matrix
    full = np.full((N, K), np.nan, dtype="float32")

    for p in parts:
        print(f"[combine] Reading {p}")
        with xr.open_dataset(p) as ds:
            nodes = ds["node"].values.astype(np.int64)
            st = ds["station"].values.astype(np.int64)
            if not np.array_equal(st, station_ids):
                raise RuntimeError(f"Station ordering mismatch in {p}")
            vals = ds["coef"].values.astype("float32")  # shape (n_block, K)
            full[nodes, :] = vals

    # Write final combined file
    out_nc = "coef_full_mhhw_matlabAligned_multi.nc"
    ds_out = xr.Dataset(
        data_vars=dict(coef=(("node", "station"), full)),
        coords=dict(
            node=("node", np.arange(N, dtype=np.int64)),
            station=("station", station_ids),
        ),
        attrs=dict(
            note="Combined per-shard outputs; MlW MATLAB-aligned; multi-station"
        ),
    )

    ds_out.to_netcdf(out_nc, encoding={"coef": {"dtype": "float32"}})
    nan_count = int(np.isnan(full).sum())
    print(f"[combine] Wrote {out_nc}")
    print(f"[combine] Combined NaNs: {nan_count} of {full.size} | shape={full.shape}")


if __name__ == "__main__":
    main()

