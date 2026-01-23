#!/usr/bin/env python3
"""
Combine all  xx_hh_ll_mp*.nc into 3 NetCDF files with fixed shapes:

  TADxx.nc  : datums(datum=7, station=N_total)
  TADhh.nc  : Htime_Hval(station=N_total, high_low=480)  # [0:239]=times, [240:479]=values
  TADll.nc  : Ltime_Lval(station=N_total, high_low=480)  # [0:239]=times, [240:479]=values

Notes
- We PRESERVE the original encoding:
  * Highs: 'higher-highs' are stored as (value + 100)
  * Lows : 'lower-lows'   are stored as (value - 100)
  * Unused slots are filled with -9.9999004e+03
  * station are mesh nodes
- This script DOES NOT decode markers or strip padding; correlation coefficient code(S3) should decode before analysis.

How to run
python S2_combine_TAD_mp.py --pattern "xx_hh_ll_mp*.nc" --out-prefix "" --chunk-stations 8192
Chunk stations is related to Pacific grid and SECOFS  8192, Alaska 16384

Created in Matlab by Liujua.Tang@noaa.gov 01/27/2020
Converted to python and changed the output to the netcdf format instead of .mat by  mojgan.rostaminia@noaa.gov, 09/08/2025
"""

#!/usr/bin/env python3

import argparse, glob, sys, re
import numpy as np
from netCDF4 import Dataset

# Exact MATLAB Fill Value
FILLV = np.float32(-9999.9004)
CAP   = 240           

def numerical_sort_key(s):
    numbers = re.findall(r'\d+', s)
    return int(numbers[-1]) if numbers else s

def scan_total_stations(file_list):
    total = 0
    widths = set()
    for fn in file_list:
        with Dataset(fn, "r") as ds:
            if "stationN" not in ds.variables:
                raise ValueError(f"{fn}: missing 'stationN'")
            total += ds.variables["stationN"].shape[0]
            wH = ds.variables["Htime_Hval"].shape[1]
            wL = ds.variables["Ltime_Lval"].shape[1]
            widths.add((wH, wL))
    for (wH, wL) in widths:
        if (wH != 2*CAP) or (wL != 2*CAP):
            raise ValueError(f"Expected width {2*CAP}; got H={wH}, L={wL}")
    return total

def create_out_TADxx(path, N_total, chunk_stations):
    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("datum", 7)
    nc.createDimension("station", None)
    v_sta = nc.createVariable("stationN", "i4", ("station",))
    
    # Write raw values directly (no _FillValue attribute)
    v_dat = nc.createVariable(
        "datums", "f4", ("datum", "station"),
        zlib=True, complevel=4,
        chunksizes=(7, min(chunk_stations, max(1, N_total)))
    )
    nc.title = "Combined tidal datums (7 x N_stations)"
    return nc, v_sta, v_dat

def create_out_hl(path, var_name, title, N_total, chunk_stations):
    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("high_low", 2*CAP)
    nc.createDimension("station", None)
    v_sta = nc.createVariable("stationN", "i4", ("station",))
    v_hl  = nc.createVariable(
        var_name, "f4", ("station", "high_low"),
        zlib=True, complevel=4, 
        chunksizes=(min(chunk_stations, max(1, N_total)), 2*CAP)
    )
    nc.title = title
    return nc, v_sta, v_hl

def safe_fill(data):
    """
    Standard fill: Replace Masked/NaN with FILLV.
    """
    if np.ma.is_masked(data):
        data = data.filled(FILLV)
    data = np.array(data, dtype=np.float32)
    data[np.isnan(data)] = FILLV
    return data

def apply_zero_patch(datums):
    """
    detects nodes where ALL 7 datums are exactly 0.0 and replaces them with FILLV.
    Input shape: (7, N_stations) or (N_stations, 7)
    """
    # Ensure shape is (7, N) for processing
    transposed = False
    if datums.shape[0] != 7:
        datums = datums.T
        transposed = True
    
    # Check for exact 0.0 across all 7 datums (axis 0)
    # returns boolean array of shape (N,)
    zero_mask = np.all(datums == 0.0, axis=0)
    
    n_zeros = np.sum(zero_mask)
    if n_zeros > 0:
        # Fill these columns with FILLV
        datums[:, zero_mask] = FILLV
    
    # Return to original shape
    if transposed:
        datums = datums.T
        
    return datums

def normalize_datums(dat):
    if dat.ndim != 2 or 7 not in dat.shape:
        raise ValueError("datums must be 2-D and have a size-7 axis")
    if dat.shape[1] == 7:
        out = dat.astype(np.float32, copy=False)
    else:
        out = dat.T.astype(np.float32, copy=False)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="xx_hh_ll_mp*.nc", help="glob for input part files")
    ap.add_argument("--out-prefix", default="", help="prefix for output files")
    ap.add_argument("--chunk-stations", type=int, default=8192, help="station chunk size")
    args = ap.parse_args()

    files = glob.glob(args.pattern)
    files.sort(key=numerical_sort_key)
    if not files:
        sys.exit(f"No files match: {args.pattern}")

    print(f"Found {len(files)} files.")
    N_total = scan_total_stations(files)
    print(f"Total stations: {N_total}")

    path_xx = f"{args.out_prefix}TADxx.nc"
    path_hh = f"{args.out_prefix}TADhh.nc"
    path_ll = f"{args.out_prefix}TADll.nc"

    nc_xx, sta_xx, dat_xx = create_out_TADxx(path_xx, N_total, args.chunk_stations)
    nc_hh, sta_hh, hl_hh  = create_out_hl(path_hh, "Htime_Hval", "Combined highs", N_total, args.chunk_stations)
    nc_ll, sta_ll, hl_ll  = create_out_hl(path_ll, "Ltime_Lval", "Combined lows", N_total, args.chunk_stations)

    try:
        offset = 0
        for fn in files:
            print(f"Processing {fn} ...")
            with Dataset(fn, "r") as ds:
                stationN = np.array(ds.variables["stationN"][:], dtype=np.int32)
                
                # 1. Read & Standard Fill (NaNs -> FILLV)
                datums = safe_fill(ds.variables["datums"][:])
                H      = safe_fill(ds.variables["Htime_Hval"][:])
                L      = safe_fill(ds.variables["Ltime_Lval"][:])
                
                # 2. Apply Zero Patch (0.0 -> FILLV)
                # Only for datums, as H/L might legitimately have 0.0 times
                datums = apply_zero_patch(datums)
                
                datums = normalize_datums(datums)

            nsta = stationN.shape[0]

            # Write
            sta_slice = slice(offset, offset + nsta)
            sta_xx[sta_slice] = stationN
            sta_hh[sta_slice] = stationN
            sta_ll[sta_slice] = stationN

            dat_xx[:, sta_slice] = datums.T
            hl_hh[sta_slice, :] = H
            hl_ll[sta_slice, :] = L

            offset += nsta
            
        print(f"Done. Wrote: {path_xx}")

    finally:
        nc_xx.close()
        nc_hh.close()
        nc_ll.close()

    # --- AUTO-VERIFICATION ---
    print("\n--- Verifying Output (Index 0, 130) ---")
    try:
        with Dataset(path_xx, 'r') as nc:
            nc.variables['datums'].set_auto_mask(False)
            val = nc.variables['datums'][0, 130]
            
            print(f"Value at [0, 130]: {val}")
            
            if np.isclose(val, FILLV, atol=1e-5):
                print("✅ SUCCESS: 0.0 was successfully replaced with -9999.9!")
            elif val == 0.0:
                 print("❌ FAIL: Value is still 0.0.")
            else:
                 print(f"❓ UNEXPECTED: {val}")

    except Exception as e:
        print(f"Verification Error: {e}")

if __name__ == "__main__":
   main()
