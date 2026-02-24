#!/usr/bin/env python3
"""
verify_phase3_full_audit.py
COMPLETE AUDIT:
1. CSV Source (RMS)
2. Phase 2 Input (uncer)
3. Phase 3 Output (d, dd)
4. Jackknife Validation (d_jn, dd_jn)
"""

import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ================= CONFIGURATION =================
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
INDATA_DIR = os.path.join(BASE_DIR, "indata")
PRE_DIR = os.path.join(BASE_DIR, "pre_process")
OUT_DIR = os.path.join(BASE_DIR, "out_nc")

# Files
CSV_FILE = os.path.join(INDATA_DIR, "Obs_Pacific_SVU.csv")
PHASE2_FILE = os.path.join(PRE_DIR, "SVU_input_mllw_Pac_SAL.nc")
PHASE3_FILE = os.path.join(OUT_DIR, "d_dd_diaPA_coef4_mllw_Pac_SAL_Pac_SAL.nc")

DATUM = "mllw"
STN_IDX = 0  # Testing the first station
# =================================================

def main():
    print(f"--- Full Phase 3 Audit (inc. Jackknife) for {DATUM} ---")
    
    # ---------------------------------------------------------
    # STEP 1: READ CSV (Source of RMS)
    # ---------------------------------------------------------
    print(f"\n[1] CSV Source")
    if not os.path.exists(CSV_FILE):
        sys.exit("   Error: CSV file not found.")
        
    df = pd.read_csv(CSV_FILE)
    df.columns = df.columns.str.lower()
    
    # Filter valid rows to find the correct Station 0
    valid_rows = []
    for i, row in df.iterrows():
        val = row[DATUM]
        r = row['rms']
        if not np.isnan(val) and val != 0 and not (np.isnan(r) or r > 999 or r < -9):
            valid_rows.append(row)
            
    target_row = valid_rows[STN_IDX]
    csv_id = int(target_row['id'])
    csv_rms = float(target_row['rms'])
    csv_val = float(target_row[DATUM])
    
    print(f"   Station ID:   {csv_id}")
    print(f"   Obs Value:    {csv_val:.4f}")
    print(f"   RMS (Uncer):  {csv_rms:.6f}")

    # ---------------------------------------------------------
    # STEP 2: PHASE 2 INPUT (Verification)
    # ---------------------------------------------------------
    print(f"\n[2] Phase 2 Input (The Transfer)")
    with Dataset(PHASE2_FILE, 'r') as nc:
        nc_uncer = nc.variables['uncer'][STN_IDX]
        nc_estn = nc.variables['estn'][STN_IDX]
        nc_ostn = nc.variables['ostn'][STN_IDX]
        node_idx = nc.variables['bcpsi'][STN_IDX] - 1
        m_error = float(nc.variables['m_error'][:])

    print(f"   NetCDF uncer: {nc_uncer:.6f}")
    print(f"   NetCDF ostn:  {nc_ostn:.4f}")
    
    if abs(csv_rms - nc_uncer) < 1e-6:
        print("   ✅ CSV RMS matched NetCDF uncer.")
    else:
        print(f"   ❌ RMS Mismatch! ({csv_rms} vs {nc_uncer})")

    # ---------------------------------------------------------
    # STEP 3: PHASE 3 OUTPUT (Final Result)
    # ---------------------------------------------------------
    print(f"\n[3] Phase 3 Output (The Result)")
    with Dataset(PHASE3_FILE, 'r') as nc:
        # Final Grid Values
        final_d = nc.variables['d'][node_idx]
        final_dd = nc.variables['dd'][node_idx] # Sigma
        
        # Jackknife Values (Station Array)
        d_jn = nc.variables['d_jn'][STN_IDX]
        dd_jn = nc.variables['dd_jn'][STN_IDX] # Sigma without this station
        
        # Change Flag (Did Jackknife reject it?)
        change_flag = nc.variables['change'][STN_IDX]

    print(f"   Final d (Datum):     {final_d:.4f}")
    print(f"   Final dd (Sigma):    {final_dd:.6f}")
    print(f"   Jackknife d_jn:      {d_jn:.4f}")
    print(f"   Jackknife dd_jn:     {dd_jn:.6f}")
    print(f"   Change Flag:         {change_flag} (1=Reweighted, 0=Normal)")

    # ---------------------------------------------------------
    # STEP 4: LOGIC CHECKS
    # ---------------------------------------------------------
    print(f"\n[4] Logic Analysis")
    
    # A. Did the final result match the observation?
    diff_final = abs(final_d - nc_ostn)
    print(f"   A. Final Error (|d - obs|): {diff_final:.6f}")
    if diff_final < nc_uncer:
        print("      ✅ Excellent. Final datum is within RMS of the observation.")
    else:
        print("      ⚠️ Final datum drifted from observation (likely due to smoothing).")

    # B. Jackknife Analysis (The "What If" Scenario)
    diff_jn = abs(d_jn - nc_ostn)
    print(f"   B. Jackknife Delta (|d_jn - obs|): {diff_jn:.6f}")
    
    if diff_jn < (nc_uncer * 3):
        print("      ✅ Neighbors support this station (Prediction matches Obs).")
    else:
        print("      ⚠️ OUTLIER WARNING: Neighbors predict a very different value!")
        print("         (This station is unique or disagrees with nearby stations.)")
        
    # C. Uncertainty Comparison
    print(f"   C. Uncertainty Hierarchy:")
    print(f"      Background: {m_error:.6f} (Least confident)")
    print(f"      Jackknife:  {dd_jn:.6f} (Prediction confidence)")
    print(f"      Final:      {final_dd:.6f} (Most confident)")
    
    if final_dd < dd_jn:
        print("      ✅ Adding this station improved confidence (Final < Jackknife).")
    else:
        print("      ❌ Something is wrong. Adding data should reduce uncertainty.")

if __name__ == "__main__":
    main()
