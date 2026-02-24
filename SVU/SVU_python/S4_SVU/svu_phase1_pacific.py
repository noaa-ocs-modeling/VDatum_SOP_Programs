#!/usr/bin/env python3
"""
svu_phase1_pacific_v4.py
Integrates Phase 0 results.
Logic:
1. Load "nanok" (Dry) and "cutok" (Stuck) from indata/Dry_Lcut_PA_SAL.nc
2. Identify "Holes" (-9999) from TADxx.nc
3. Identify "Spikes" (> 0.001) from TADxx.nc
4. Combine ALL Bad Nodes -> Interpolate.
"""

import os
import sys
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import NearestNDInterpolator

# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
FORT14_FILE = os.path.join(BASE_DIR, "indata/fort.14")
TAD_FILE    = "/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_08182025_TAD/TADxx.nc"
FORT63_MASK = os.path.join(BASE_DIR, "indata/Dry_Lcut_PA_SAL.nc") # <--- Output of Phase 0
OUT_DIR     = os.path.join(BASE_DIR, "pre_process")

DATUMS = ["mhhw", "mhw", "msl", "mtl", "dtl", "mlw", "mllw"]
DATUM_INDICES = {"mhhw": 0, "mhw": 1, "msl": 2, "mtl": 3, "dtl": 4, "mlw": 5, "mllw": 6}
SPIKE_THRESHOLD = 0.001 

# ============================================================
# UTILITIES
# ============================================================
def read_fort14(f14_path):
    print(f"Reading coordinates from {f14_path}...")
    with open(f14_path, 'r') as f:
        f.readline()
        line2 = f.readline().split()
        NP = int(line2[1])
    import pandas as pd
    df = pd.read_csv(f14_path, delim_whitespace=True, header=None, 
                     skiprows=2, nrows=NP, usecols=[1, 2], names=['lon', 'lat'])
    return df['lon'].values, df['lat'].values

def save_nc(path, data, name="z"):
    if os.path.exists(path): os.remove(path)
    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("node", len(data))
        v = nc.createVariable(name, "f4", ("node",))
        v[:] = data

def save_mask_nc(path, mask_array):
    if os.path.exists(path): os.remove(path)
    with Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("node", len(mask_array))
        v = nc.createVariable("node", "i4", ("node",))
        # Save as 1=Valid, 0=Invalid (Interpolated)
        v[:] = mask_array.astype(np.int32)

def clean_datum(datum_name, raw_data, lon, lat, out_dir, fort63_bad_mask):
    print(f"\n--- Cleaning {datum_name} ---")
    data = raw_data.copy()
    data[data < -1000] = np.nan # Mark original holes
    
    # 1. SPIKE FILTER (For Low Water Datums)
    if datum_name in ["mlw", "mllw"]:
        spike_mask = data > SPIKE_THRESHOLD
        n_spikes = np.sum(spike_mask)
        if n_spikes > 0:
            print(f"   Detected {n_spikes} positive spikes. Marking as NaN.")
            data[spike_mask] = np.nan

    # 2. FORT.63 FILTER (Dry/Stuck Nodes)
    if fort63_bad_mask is not None:
        # fort63_bad_mask is True where node is BAD
        bad_idx = np.where(fort63_bad_mask)[0]
        
        # Check how many NEW bad nodes this adds
        existing_nans = np.isnan(data)
        new_bad = np.sum(~existing_nans[bad_idx])
        
        if new_bad > 0:
            print(f"   Merging Phase 0 defects: Removing {new_bad} additional 'stuck/dry' nodes.")
            data[bad_idx] = np.nan

    # 3. INTERPOLATE
    is_valid = ~np.isnan(data)
    n_nans = np.sum(~is_valid)
    
    # Save Mask
    mask_path = os.path.join(out_dir, f"mdatum_nanloc_{datum_name}.nc")
    save_mask_nc(mask_path, is_valid)
    print(f"   Total Invalid Count: {n_nans}")
    print(f"   Saved Mask: {mask_path}")

    if n_nans > 0:
        print(f"   Interpolating {n_nans} nodes...")
        interp = NearestNDInterpolator(list(zip(lon[is_valid], lat[is_valid])), data[is_valid])
        filled_values = interp(list(zip(lon[~is_valid], lat[~is_valid])))
        data[~is_valid] = filled_values
    
    # 4. SAVE CLEAN FILE
    out_path = os.path.join(out_dir, f"dm_{datum_name}.nc")
    save_nc(out_path, data, "z")

# ============================================================
# MAIN
# ============================================================
def main():
    if not os.path.exists(OUT_DIR): os.makedirs(OUT_DIR)
    lon, lat = read_fort14(FORT14_FILE)

    # 1. Load Phase 0 Mask
    fort63_bad_mask = None
    if os.path.exists(FORT63_MASK):
        print(f"Loading Phase 0 defects from {FORT63_MASK}...")
        with Dataset(FORT63_MASK, 'r') as nc:
            # Load both masks (1=Bad)
            nanok = nc.variables['nanok'][:]
            cutok = nc.variables['cutok'][:]
            
            # Combine: True if EITHER nanok OR cutok is 1
            fort63_bad_mask = (nanok == 1) | (cutok == 1)
            
            print(f"   Found {np.sum(fort63_bad_mask)} bad nodes (Dry + Stuck).")
    else:
        print(f"WARNING: {FORT63_MASK} not found. Did you run svu_phase0_pacific.py?")

    # 2. Process Datums
    print(f"Opening {TAD_FILE}...")
    with Dataset(TAD_FILE, 'r') as nc:
        datums_arr = nc.variables['datums'][:]
        if datums_arr.shape[0] != 7: datums_arr = datums_arr.T
        
        for dname in DATUMS:
            clean_datum(dname, datums_arr[DATUM_INDICES[dname], :len(lon)], 
                        lon, lat, OUT_DIR, fort63_bad_mask)

    print("\n=== Phase 1 v4 Complete ===")

if __name__ == "__main__":
    main()
