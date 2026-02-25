#!/usr/bin/env python
# coding: utf-8

# There are two scripts at "C:\Users\mojgan.rostaminia\Documents\Hawaii_Pacific\IOC_data\compute_uncertainy_IOC.ipynb" and 
# "C:\Users\mojgan.rostaminia\Documents\Hawaii_Pacific\DART_data\compute_uncertainy_DART.ipynb"  to do this sepratly. The result are identical

# In[7]:


import os
import numpy as np
import pandas as pd
from pyproj import Geod

# =========================
# CONFIGURATION & CONSTANTS
# =========================
K_HOURS = 6.21      # 12.42/2 (timing proxy scale)
M2FT = 3.28084      # Meters to Feet conversion
_GEOD = Geod(ellps="WGS84")

# =========================
# HELPER FUNCTIONS
# =========================
def norm_id(x):
    """Normalizes station IDs by stripping whitespace and trailing decimals."""
    return str(x).strip().replace(".0", "")

def _dist_km_karney_pyproj(lon1, lat1, lon2, lat2):
    """Calculates ellipsoidal (WGS-84) distance. Vectorized. Returns km."""
    _, _, s12 = _GEOD.inv(lon1, lat1, lon2, lat2)
    return s12 / 1000.0

# =========================
# MAIN PIPELINE FUNCTION
# =========================
def compute_uncertainty(xlsx_path, station_csv, out_csv, method='3mo'):
    """
    method: '3mo' for Bodnar 3-month (Pacific), 'long' for Bodnar long-record (Pacific >= 12mo)
    """
    print(f"\n--- Processing: {os.path.basename(xlsx_path)} (Method: {method}) ---")
    
    # 1. READ DATUMS
    ctl = pd.read_excel(xlsx_path, sheet_name="ControlStation")
    ctl.columns = [str(c).strip() for c in ctl.columns]
    
    sub = pd.read_excel(xlsx_path, sheet_name="SubordinateStation")
    sub.columns = [str(c).strip() for c in sub.columns]
    var_col = sub.columns[0]
    sub[var_col] = sub[var_col].astype(str).str.strip()
    sub = sub.set_index(var_col)
    sub.columns = [norm_id(c) for c in sub.columns]
    
    # Normalize IDs & Convert Numerics
    for c in ["Station ID_ctl", "Station ID_sub"]:
        if c in ctl.columns:
            ctl[c] = ctl[c].apply(norm_id)
            
    for c in ["MN", "DHQ", "DLQ"]:
        if c in ctl.columns:
            ctl[c] = pd.to_numeric(ctl[c], errors="coerce")
            
    ctl = ctl.dropna(subset=["Station ID_ctl", "Station ID_sub", "MN", "DHQ", "DLQ"]).copy()
    
    # Base pairs dataframe
    pairs = ctl[["Station ID_ctl", "Station ID_sub"]].copy()
    
    # 2. FETCH MEAN RANGES (MN)
    MN_sub_series = pd.to_numeric(sub.loc["MN"], errors="coerce")
    MN_sub_lookup = MN_sub_series.to_dict()
    pairs["MN_sub_m"] = pairs["Station ID_sub"].map(MN_sub_lookup)
    
    MN_ctl_lookup = ctl.groupby("Station ID_ctl")["MN"].median().to_dict()
    pairs["MN_ctl_m"] = pairs["Station ID_ctl"].map(MN_ctl_lookup)
    
    # Convert to Feet for Bodnar Equations
    pairs["MN_sub_ft"] = pairs["MN_sub_m"] * M2FT
    pairs["MN_ctl_ft"] = pairs["MN_ctl_m"] * M2FT
    
    # Mean Range Ratio (MNR)
    pairs["MNR"] = (pairs["MN_sub_ft"] - pairs["MN_ctl_ft"]).abs() / pairs["MN_ctl_ft"]
    
    # Standard columns for both outputs
    out_cols = [
        "Station ID_ctl", "Station ID_sub", 
        "MN_sub_m", "MN_ctl_m", "MN_sub_ft", "MN_ctl_ft", "MNR"
    ]

    # -----------------------------------------
    # CONDITIONAL LOGIC BASED ON METHOD
    # -----------------------------------------
    if method == '3mo':
        # COMPUTE TIMING PROXY (T)
        ctl["r"] = (ctl["DHQ"] - ctl["DLQ"]) / ctl["MN"]
        pairs["T_hr"] = (K_HOURS * ctl["r"]).abs()
        
        # COMPUTE ELLIPSOIDAL DISTANCE
        meta = pd.read_csv(station_csv)
        meta.columns = [str(c).strip() for c in meta.columns]
        meta["STATION_ID"] = meta["STATION_ID"].apply(norm_id)
        meta["ST_LATITUDE"]  = pd.to_numeric(meta["ST_LATITUDE"], errors="coerce")
        meta["ST_LONGITUDE"] = pd.to_numeric(meta["ST_LONGITUDE"], errors="coerce")
        meta = meta.dropna(subset=["STATION_ID", "ST_LATITUDE", "ST_LONGITUDE"])
        
        lat_lookup = dict(zip(meta["STATION_ID"], meta["ST_LATITUDE"]))
        lon_lookup = dict(zip(meta["STATION_ID"], meta["ST_LONGITUDE"]))
        
        pairs["lat_ctl"] = pairs["Station ID_ctl"].map(lat_lookup)
        pairs["lon_ctl"] = pairs["Station ID_ctl"].map(lon_lookup)
        pairs["lat_sub"] = pairs["Station ID_sub"].map(lat_lookup)
        pairs["lon_sub"] = pairs["Station ID_sub"].map(lon_lookup)
        
        pairs_ok = pairs.dropna(subset=["lat_ctl", "lon_ctl", "lat_sub", "lon_sub"]).copy()
        
        pairs_ok["dist_km"] = _dist_km_karney_pyproj(
            pairs_ok["lon_ctl"].values, pairs_ok["lat_ctl"].values,
            pairs_ok["lon_sub"].values, pairs_ok["lat_sub"].values
        )
        pairs_ok["dist_nm"] = pairs_ok["dist_km"] / 1.852
        pairs_ok["D_sqrt_nm"] = np.sqrt(pairs_ok["dist_nm"])
        
        # Merge distances back to main dataframe
        pairs = pairs.merge(
            pairs_ok[["Station ID_ctl", "Station ID_sub", "dist_km", "dist_nm", "D_sqrt_nm"]],
            on=["Station ID_ctl", "Station ID_sub"],
            how="left"
        )
        
        # Bodnar 3-Month Model
        pairs["sigma_ft_3mo"] = (0.0043 * pairs["T_hr"] + 0.0036 * pairs["D_sqrt_nm"] + 0.0255 * pairs["MNR"] + 0.029)
        pairs["sigma_m_3mo"] = pairs["sigma_ft_3mo"] * 0.3048
        
        # Reorder columns to put Time and Distance up front
        out_cols = ["Station ID_ctl", "Station ID_sub", "T_hr", "dist_km", "dist_nm", "D_sqrt_nm"] + out_cols[2:] + ["sigma_ft_3mo", "sigma_m_3mo"]
        
    elif method == 'long':
        # Bodnar Long-Record Model (>= 12 months)
        pairs["SRSMN"] = np.sqrt(pairs["MN_sub_ft"] + pairs["MN_ctl_ft"])
        pairs["sigma_ft_long"] = (0.0045 * pairs["SRSMN"] + 0.0128 * pairs["MNR"] + 0.025)
        pairs["sigma_m_long"] = pairs["sigma_ft_long"] * 0.3048
        
        # Append long-record specific columns
        out_cols += ["SRSMN", "sigma_ft_long", "sigma_m_long"]
        
    else:
        raise ValueError("Error: Method must be '3mo' or 'long'.")
    
    # 3. SAVE SPECIFIC OUTPUT
    final_df = pairs[out_cols]
    final_df.to_csv(out_csv, index=False)
    print(f"Saved results to: {out_csv}")
    
    return final_df

# =========================
# EXECUTION
# =========================
if __name__ == "__main__":
    # Common Station Coordinates
    STATION_CSV = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/stations/Pacific_station_final/Obs_Pacific.csv"
    
    # 1. Process DART Data (Using 3-Month Method)
    dart_xlsx = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/DART_data/DartBuoy_Datums.xlsx"
    dart_out = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/DART_data/Bodnar_uncertainties_DART.csv"
    
    # Uncomment to run:
    compute_uncertainty(dart_xlsx, STATION_CSV, dart_out, method='3mo')
    
    # 2. Process IOC Data (Using Long-Record Method)
    ioc_xlsx = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/IOC_data/IOC_Datums.xlsx"
    ioc_out = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/IOC_data/Bodnar_uncertainties_IOC.csv"
    
    # Uncomment to run:
    compute_uncertainty(ioc_xlsx, STATION_CSV, ioc_out, method='long')


# In[ ]:




