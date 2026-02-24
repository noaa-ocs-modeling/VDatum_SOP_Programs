#!/usr/bin/env python3
import os
import argparse
import numpy as np
from netCDF4 import Dataset
import sys

# ------------------------------------------------------------------------------
# 1. CORE MATH FUNCTIONS
# ------------------------------------------------------------------------------

def get_station_residuals(estn, uncer, dm, W01, covar0, bcpsi_idx, idatum_idx):
    """
    Calculates residuals at stations. Returns weighted matrices for d/dd calculation.
    """
    estn = estn.reshape(-1, 1)
    uncer = uncer.reshape(-1, 1)
    dm = dm.reshape(-1, 1)
    
    HPH = covar0[bcpsi_idx, :] 
    m = HPH.shape[0]

    # Construct R matrix (Weighted with W01)
    # R is diagonal of squared uncertainty
    R_diag = np.zeros(len(uncer))
    for i in range(len(uncer)):
        R_diag[i] = W01[i, i] * (uncer[i, 0]**2)
    R = np.diag(R_diag)
    
    RHPH = R + HPH
    
    ones_col = np.ones((m, 1))
    RHPH1 = np.block([
        [RHPH, ones_col],
        [ones_col.T, 0.0]
    ])
    
    try:
        invRHPH1 = np.linalg.inv(RHPH1)
    except np.linalg.LinAlgError:
        invRHPH1 = np.linalg.pinv(RHPH1)

    HPH1 = np.hstack([HPH, ones_col])
    G_stations_full = HPH1 @ invRHPH1
    G_stations = G_stations_full[:, :m]
    
    dd_stations = G_stations @ estn
    
    dm_stations = dm[bcpsi_idx]
    d_stations = dm_stations + dd_stations
    
    # Sign Correction Logic
    if idatum_idx < 4:
        with np.errstate(invalid='ignore'): 
            sign_change_mask = (np.sign(d_stations) != np.sign(dm_stations)) & (dm_stations != 0) & (~np.isnan(dm_stations))
        
        loc_sign = np.where(sign_change_mask)[0]
        if len(loc_sign) > 0:
            d_stations[loc_sign] = dm_stations[loc_sign]
            dd_stations[loc_sign] = 0.0
            
    return d_stations, invRHPH1, R, HPH

def calculate_full_grid(covar0, invRHPH1, estn, dm, m, uncer, HPH, m_error, idatum_idx):
    """
    Calculates full grid d, dd, and uncertainty (diaPA).
    FIX: Uses unweighted R for diaPA calculation to match MATLAB.
    """
    print("Calculating full grid...", flush=True)
    n = covar0.shape[0]
    
    # --- 1. Calculate Weights (G) and Corrections (dd) ---
    inv_top = invRHPH1[:m, :] 
    G_full = covar0 @ inv_top
    
    last_row_inv = invRHPH1[m, :] 
    G_full += last_row_inv 
    
    G = G_full[:, :m]
    
    dd = G @ estn
    d = dm + dd
    
    del G_full 
    
    # --- 2. Sign Correction (Full Grid) ---
    if idatum_idx < 4:
        with np.errstate(invalid='ignore'):
            sign_change_mask = (np.sign(d) != np.sign(dm)) & (dm != 0) & (~np.isnan(dm))
        loc_sign = np.where(sign_change_mask)[0]
        if len(loc_sign) > 0:
            print(f"Applying sign correction to {len(loc_sign)} grid points.", flush=True)
            d[loc_sign] = dm[loc_sign]
            dd[loc_sign] = 0.0

    # --- 3. Calculate Uncertainty (diaPA) ---
    print(f"Calculating diaPA (using runtime m_error={m_error:.6f})...", flush=True)
    
    diaP = (m_error**2) * np.ones((n, 1))
    
    term2 = np.sum(G * covar0, axis=1, keepdims=True)
    
    # Reconstruct R WITHOUT weights for uncertainty calculation
    R_pure = np.diag((uncer.flatten()**2))
    
    temp_small = R_pure + HPH
    G_temp = G @ temp_small
    term3 = np.sum(G_temp * G, axis=1, keepdims=True)
    
    del G_temp
    del G 
    
    variance_PA = diaP - 2.0 * term2 + term3
    
    # Handle Negative Variance (Complex Numbers)
    diaPA = np.zeros(variance_PA.shape, dtype=np.complex128)
    
    mask_pos = variance_PA >= 0
    diaPA[mask_pos] = np.sqrt(variance_PA[mask_pos])
    
    mask_neg = ~mask_pos
    if np.any(mask_neg):
        print(f"Notice: {np.sum(mask_neg)} points have negative variance (stored in diaPA_imag).", flush=True)
        diaPA[mask_neg] = np.sqrt(np.abs(variance_PA[mask_neg])) * 1j

    return d, dd, diaPA

# ------------------------------------------------------------------------------
# 2. FILE PROCESSING
# ------------------------------------------------------------------------------

def process_datum(runid, datum, path_pre, path_out, array_idx, coef_method=4):
    print(f'----SVU for {datum}-------------', flush=True)
    
    # Input File: NetCDF from Phase 2
    file_input = os.path.join(path_pre, f"SVU_input_{datum}_{runid}.nc")
    
    if not os.path.exists(file_input):
        print(f"CRITICAL: Input NetCDF file not found: {file_input}", flush=True)
        return

    print(f"Loading input from NetCDF: {file_input}...", flush=True)
    
    try:
        with Dataset(file_input, 'r') as nc:
            nc.set_auto_mask(False)
            estn   = nc.variables['estn'][:].astype(np.float64)
            uncer  = nc.variables['uncer'][:].astype(np.float64)
            bcpsi  = nc.variables['bcpsi'][:].astype(int)
            ostn   = nc.variables['ostn'][:].astype(np.float64)
            dm     = nc.variables['dm'][:].astype(np.float64)
            covar0 = nc.variables['covar0'][:].astype(np.float64)
            
    except Exception as e:
        print(f"Error reading NetCDF file: {e}")
        return

    # Reshape for math functions
    estn   = estn.reshape(-1, 1)
    uncer  = uncer.reshape(-1, 1)
    dm     = dm.reshape(-1, 1)
    
    if ostn is not None:
        ostn = ostn.flatten()
    else:
        ostn = np.zeros_like(uncer).flatten()

    print(f"Data Shapes -> dm: {dm.shape}, covar: {covar0.shape}", flush=True)

    # Convert 1-based indices to 0-based
    bcpsi_idx = (bcpsi.flatten() - 1).astype(int)
    
    num_stations = len(uncer)
    W01 = np.eye(num_stations)
    
    # Jackknife Loop
    for mc in range(50):
        d_stations, invRHPH1, R, HPH = get_station_residuals(
            estn, uncer, dm, W01, covar0, bcpsi_idx, array_idx
        )
        change = np.zeros(num_stations, dtype=int)
        d_flat = d_stations.flatten()
        
        for i in range(num_stations):
            if np.isnan(d_flat[i]) or np.isnan(ostn[i]):
                continue
            diff = abs(d_flat[i] - ostn[i])
            curr_uncer = uncer[i].item()
            threshold = max(min(0.01, curr_uncer), 1e-6)
            if diff > threshold:
                W01[i, i] = W01[i, i] / np.sqrt(2.0)
                change[i] = 1
        
        if np.sum(change) == 0:
            break

    m_error = np.std(estn, ddof=1)
    print(f"Recalculated m_error: {m_error}", flush=True)

    m = HPH.shape[0]
    
    # Calculate Full Grid
    d, dd, diaPA = calculate_full_grid(covar0, invRHPH1, estn, dm, m, uncer, HPH, m_error, array_idx)
    
    # Extract Jackknife results for stations
    d_jn = d[bcpsi_idx]
    dd_jn = dd[bcpsi_idx]
    diaPA_jn = diaPA[bcpsi_idx]

    # --- OUTPUT: Save as NetCDF ---
    outfile = os.path.join(path_out, f'd_dd_diaPA_{datum}_{runid}.nc')
    print(f"Saving to {outfile}", flush=True)

    try:
        with Dataset(outfile, "w", format="NETCDF4") as nc:
            # Dimensions
            nc.createDimension("node", len(dm))
            nc.createDimension("station", len(change))
            nc.createDimension("one", 1)

            # Variables - All f8 (double)
            # Full Grid
            nc.createVariable("d", "f8", ("node",))[:] = d.flatten()
            nc.createVariable("dd", "f8", ("node",))[:] = dd.flatten()
            nc.createVariable("dm", "f8", ("node",))[:] = dm.flatten()
            
            # Complex Variables (Split Real/Imag)
            nc.createVariable("diaPA_real", "f8", ("node",))[:] = np.real(diaPA).flatten()
            nc.createVariable("diaPA_imag", "f8", ("node",))[:] = np.imag(diaPA).flatten()

            # Station Variables
            nc.createVariable("change", "i4", ("station",))[:] = change.flatten()
            nc.createVariable("d_jn", "f8", ("station",))[:] = d_jn.flatten()
            nc.createVariable("dd_jn", "f8", ("station",))[:] = dd_jn.flatten()
            nc.createVariable("diaPA_jn_real", "f8", ("station",))[:] = np.real(diaPA_jn).flatten()
            nc.createVariable("diaPA_jn_imag", "f8", ("station",))[:] = np.imag(diaPA_jn).flatten()

            # Scalar
            nc.createVariable("final_m_error", "f8")[:] = m_error

        print(f"SUCCESS: Finished {datum}. File saved.", flush=True)
        
    except Exception as e:
        print(f"Error saving output NetCDF: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--runid', type=str, required=True)
    parser.add_argument('--path_pre', type=str, required=True)
    parser.add_argument('--path_out', type=str, required=True)
    parser.add_argument('--array_index', type=int, required=True)
    args = parser.parse_args()
    
    # Order matches SLURM array
    datums = ["mllw", "mhhw", "mhw", "mlw", "mtl", "dtl"] 
    
    if 0 <= args.array_index < len(datums):
        process_datum(args.runid, datums[args.array_index], args.path_pre, args.path_out, args.array_index)
