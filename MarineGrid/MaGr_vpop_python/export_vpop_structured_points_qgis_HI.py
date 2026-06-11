#!/usr/bin/env python3
"""
export_vpop_qgis_fringe.py

Export only the coastal-fringe structured VPOP points for QGIS,
following the same method used in your marine-grid script.

Method:
- read marine grid field + box
- reconstruct lon/lat from box = [alat0, alon0, dely, delx]
- wrap longitudes > 180 to -180..180 for QGIS
- create coastal fringe mask using binary dilation
- read VPOP datum
- align VPOP shape to marine grid
- export masked points to CSV
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.ndimage import binary_dilation


def masked_to_nan(a):
    if np.ma.isMaskedArray(a):
        return np.asarray(a.filled(np.nan), dtype=np.float64)
    return np.asarray(a, dtype=np.float64)


def read_marine_grid(marine_nc: str | Path):
    marine_nc = Path(marine_nc)
    print(f"Reading marine grid: {marine_nc}")
    with Dataset(marine_nc, "r") as ds:
        field = np.asarray(ds.variables["field"][:])
        alat0, alon0, dely, delx = np.asarray(ds.variables["box"][:], dtype=np.float64)

    if field.ndim != 2:
        raise ValueError(f"{marine_nc}: expected 2D field, got {field.shape}")

    jmax, imax = field.shape
    clon_raw = alon0 + delx * np.arange(imax, dtype=np.float64)
    clat = alat0 + dely * np.arange(jmax, dtype=np.float64)

    # QGIS-friendly longitude wrap
    clon = np.where(clon_raw > 180.0, clon_raw - 360.0, clon_raw)

    return field, clon, clat


def read_vpop_field(vpop_nc: str | Path):
    vpop_nc = Path(vpop_nc)
    print(f"Reading VPOP field: {vpop_nc}")

    with Dataset(vpop_nc, "r") as ds:
        if "datum" in ds.variables:
            arr = masked_to_nan(ds.variables["datum"][:])
            if arr.ndim == 2:
                return arr, "datum"

        for vname, var in ds.variables.items():
            arr = masked_to_nan(var[:])
            if arr.ndim == 2:
                return arr, vname

    raise KeyError(f"{vpop_nc}: no usable 2D variable found")


def align_vpop_to_marine(vpop_arr: np.ndarray, marine_shape: tuple[int, int]):
    if vpop_arr.shape == marine_shape:
        return vpop_arr, "as_is"
    if vpop_arr.T.shape == marine_shape:
        return vpop_arr.T, "transposed"
    raise ValueError(
        f"VPOP shape {vpop_arr.shape} does not match marine shape {marine_shape}, "
        f"and transpose {vpop_arr.T.shape} does not match either"
    )


def build_fringe_mask(field: np.ndarray, iterations: int):
    """
    Match your marine-grid visualization method:
    - expansion layers are field >= 2
    - keep water field == 1 near expansion layers
    """
    expansion_mask = (field >= 2)
    near_shore_water_area = binary_dilation(expansion_mask, iterations=iterations)
    water_fringe_mask = (field == 1) & near_shore_water_area
    final_mask = expansion_mask | water_fringe_mask
    return final_mask


def main():
    ap = argparse.ArgumentParser(description="Export coastal-fringe structured VPOP points for QGIS")
    ap.add_argument("--marine-nc", required=True, help="Marine grid nc")
    ap.add_argument("--vpop-nc", required=True, help="VPOP datum nc")
    ap.add_argument("--out-csv", required=True, help="Output CSV")
    ap.add_argument("--iterations", type=int, default=10,
                    help="binary_dilation iterations for near-shore water fringe")
    ap.add_argument("--finite-only", action="store_true",
                    help="Keep only finite vpop values")
    args = ap.parse_args()

    marine_field, clon, clat = read_marine_grid(args.marine_nc)
    vpop_raw, vname = read_vpop_field(args.vpop_nc)
    vpop_field, align_mode = align_vpop_to_marine(vpop_raw, marine_field.shape)

    print("Applying coastal fringe filter...")
    final_mask = build_fringe_mask(marine_field, args.iterations)

    if args.finite_only:
        final_mask = final_mask & np.isfinite(vpop_field)

    j_idx, i_idx = np.where(final_mask)

    if len(j_idx) == 0:
        print("Error: No points found after applying coastal fringe mask.")
        return

    print(f"Extracting {len(j_idx)} points...")

    df = pd.DataFrame({
        "longitude": clon[i_idx],
        "latitude": clat[j_idx],
        "marine_value": marine_field[j_idx, i_idx],
        "vpop_value": vpop_field[j_idx, i_idx],
    })

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    print(f"\nSUCCESS: Created {out_csv}")
    print(f"Marine shape     : {marine_field.shape}")
    print(f"VPOP variable    : {vname}")
    print(f"VPOP raw shape   : {vpop_raw.shape}")
    print(f"Alignment mode   : {align_mode}")
    print(f"Longitude range  : {df['longitude'].min()} to {df['longitude'].max()}")
    print(f"Latitude range   : {df['latitude'].min()} to {df['latitude'].max()}")
    print("Action: Drag this CSV into QGIS. Set Geometry CRS to EPSG:4326.")
    print("Style by graduated colors using: vpop_value")


if __name__ == "__main__":
    main()
