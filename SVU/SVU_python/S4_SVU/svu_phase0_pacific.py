#!/usr/bin/env python3
"""
scan_fort63_defects.py
Replicates MATLAB logic: a0_check_63_nan_low_water_cut.m
Scans fort.63.nc for:
1. Dry Nodes (nanok): Value is -99999
2. Stuck Nodes (cutok): Value is constant for 3 consecutive steps (range < 1e-8)

Output: indata/Dry_Lcut_PA_SAL.nc
"""

import os
import sys
import numpy as np
from netCDF4 import Dataset

# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
FORT63_FILE = "/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/PA_6sec_SAL_08182025/fort.63.nc"
OUTPUT_FILE = os.path.join(BASE_DIR, "indata/Dry_Lcut_PA_SAL.nc")

def main():
    if not os.path.exists(FORT63_FILE):
        sys.exit(f"Error: fort.63.nc not found at {FORT63_FILE}")
        
    print(f"Scanning {FORT63_FILE}...")
    
    with Dataset(FORT63_FILE, 'r') as nc:
        # Get dimensions
        n_nodes = nc.dimensions['node'].size
        n_times = nc.dimensions['time'].size
        var_zeta = nc.variables['zeta']
        
        print(f"Nodes: {n_nodes}, Time Steps: {n_times}")
        
        # Initialize Masks (0=Good, 1=Bad)
        nanok = np.zeros(n_nodes, dtype=np.int8) # Dry
        cutok = np.zeros(n_nodes, dtype=np.int8) # Stuck/Cut
        
        # Iterate through time in chunks of 3 (Replicates MATLAB logic)
        count = 0
        
        for t in range(0, n_times, 3):
            # Ensure we have a full window of 3
            if t + 3 > n_times:
                break
                
            # Read 3 time steps for ALL nodes
            data_chunk = var_zeta[t:t+3, :] 
            
            # 1. Check for DRY (-99999)
            # Find any node that is dry in this chunk
            chunk_dry = (data_chunk < -90000)
            
            # Update global dry mask
            nanok[np.any(chunk_dry, axis=0)] = 1
            
            # 2. Check for STUCK (Low Water Cut)
            # Create working copy to handle NaNs/Dry values
            working_data = data_chunk.copy()
            working_data[chunk_dry] = np.nan
            
            # Calculate Range (Max - Min) along time axis
            dmax = np.nanmax(working_data, axis=0)
            dmin = np.nanmin(working_data, axis=0)
            diff = dmax - dmin
            
            # If diff is tiny (approx 0) and NOT all NaNs, it is stuck
            valid_window = ~np.all(np.isnan(working_data), axis=0)
            stuck_now = (np.abs(diff) < 1e-8) & valid_window
            
            cutok[stuck_now] = 1
            
            count += 3
            if count % 300 == 0:
                print(f"   Processed {count}/{n_times} steps...")

    # Summary
    n_nanok = np.sum(nanok)
    n_cutok = np.sum(cutok)
    
    print("-" * 50)
    print(f"Scan Complete.")
    print(f"Dry Nodes (nanok):    {n_nanok}")
    print(f"Stuck Nodes (cutok):  {n_cutok}")
    print("-" * 50)
    
    # Save Result
    print(f"Saving output to {OUTPUT_FILE}...")
    if os.path.exists(OUTPUT_FILE): os.remove(OUTPUT_FILE)
    
    with Dataset(OUTPUT_FILE, 'w', format="NETCDF4") as nc_out:
        nc_out.createDimension("node", n_nodes)
        
        v_nanok = nc_out.createVariable("nanok", "i1", ("node",))
        v_cutok = nc_out.createVariable("cutok", "i1", ("node",))
        
        v_nanok[:] = nanok
        v_cutok[:] = cutok

    print("Done.")

if __name__ == "__main__":
    main()
