#!/usr/bin/env python3
"""
svu_visualize_covariance_components.py
DEEP DIVE VISUALIZATION (JET + WHITE TEXT):
Reconstructs and maps the components of the VDatum Covariance Kernel.

Updates:
- Uses 'jet' colormap.
- Forces Colorbar Text to WHITE (SVG fill).
"""

import os
import sys
import numpy as np
import pandas as pd
import folium
import branca.colormap as cm
from netCDF4 import Dataset

# ================= CONFIGURATION =================
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
INDATA_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/indata/"

# 1. Final Covariance File
SVU_INPUT = os.path.join(BASE_DIR, "pre_process/SVU_input_mllw_Pac_SAL.nc")

# 2. Model Datum File
DM_FILE = os.path.join(BASE_DIR, "pre_process/dm_mllw.nc")

# 3. Observation CSV
OBS_FILE = os.path.join(INDATA_DIR, "Obs_Pacific_SVU.csv")

# 4. Component Source Directories
WDIST_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S1_Wdist/"
KCOEF_FILE = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S3_CC/Pacific/coef_full_mllw_matlabAligned_multi.nc"

STATION_IDX = 0 
THRESHOLD_PERCENT = 0.01 
FORT14_FILE = os.path.join(INDATA_DIR, "fort.14")
# =================================================

def read_fort14_coords(f14_path):
    print(f"Reading mesh coordinates from {os.path.basename(f14_path)}...")
    with open(f14_path, 'r') as f:
        f.readline()
        line2 = f.readline().split()
        NP = int(line2[1])
    
    df = pd.read_csv(f14_path, delim_whitespace=True, header=None, 
                     skiprows=2, nrows=NP, usecols=[0, 1, 2], names=['id', 'lon', 'lat'])
    return df

def generate_map(df_bubble, value_col, stn_lat, stn_lon, filename, caption):
    print(f"   Generating {filename}...")
    m = folium.Map(location=[stn_lat, stn_lon], zoom_start=9, 
                   tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                   attr='Esri World Imagery')

    min_v, max_v = df_bubble[value_col].min(), df_bubble[value_col].max()
    if max_v == min_v: max_v += 0.0001
    
    # JET COLORMAP
    colors = ['#000080', '#0000FF', '#0040FF', '#0080FF', '#00BFFF', '#00FFFF', 
              '#40FFBF', '#80FF80', '#BFFF40', '#FFFF00', '#FFBF00', '#FF8000', 
              '#FF4000', '#FF0000', '#800000']
    
    cmap = cm.LinearColormap(colors=colors, vmin=min_v, vmax=max_v, caption=caption)
    
    # CSS INJECTION: Force White Text on Legend & Ticks
    m.get_root().html.add_child(folium.Element("""
        <style>
            /* Title of the legend */
            .caption { 
                color: white !important; 
                font-weight: bold; 
                text-shadow: 1px 1px 2px black; 
                font-size: 14px;
            }
            /* Tick labels on the colorbar (SVG text) */
            .legend text { 
                fill: white !important; 
                font-weight: bold; 
                text-shadow: 1px 1px 2px black;
            }
            /* Attribution text at bottom right */
            .leaflet-control-attribution { color: #ccc !important; }
        </style>
    """))

    fg = folium.FeatureGroup(name=caption)
    
    step = 1
    if len(df_bubble) > 10000: step = int(len(df_bubble) / 5000)

    for idx in range(0, len(df_bubble), step):
        row = df_bubble.iloc[idx]
        val = row[value_col]
        
        folium.CircleMarker(
            location=[row['lat'], row['lon']],
            radius=3,
            color=cmap(val),
            fill=True,
            fill_color=cmap(val),
            fill_opacity=0.8,
            popup=f"Node: {int(row['id'])}<br>{caption}: {val:.6f}"
        ).add_to(fg)
    
    fg.add_to(m)
    folium.Marker([stn_lat, stn_lon], popup="STATION", icon=folium.Icon(color='white', icon='star', icon_color='red')).add_to(m)
    m.add_child(cmap)
    m.save(filename)

def main():
    if not os.path.exists(SVU_INPUT): sys.exit(f"Error: {SVU_INPUT} not found.")

    # 1. LOAD MAIN INFO
    print(f"Loading data for Station Index {STATION_IDX}...")
    with Dataset(SVU_INPUT, 'r') as nc:
        cov_col = nc.variables['covar0'][:, STATION_IDX]
        stn_node_idx = int(nc.variables['bcpsi'][STATION_IDX]) - 1 
        m_error = float(nc.variables['m_error'][:])
        Lr = float(nc.variables['Lr'][:])
    
    m_error_sq = m_error ** 2
    print(f"   m_error^2 (variance): {m_error_sq:.6f}")

    # 2. IDENTIFY REAL STATION ID
    df_obs = pd.read_csv(OBS_FILE)
    df_obs.columns = df_obs.columns.str.lower()
    target_row = df_obs[df_obs['node'] == (stn_node_idx + 1)]
    if len(target_row) == 0: sys.exit(f"Error: Node {stn_node_idx+1} not found in CSV.")
    real_station_id = int(target_row.iloc[0]['id'])
    print(f"   Station ID: {real_station_id}")

    # 3. FILTER BUBBLE
    max_cov = np.max(cov_col)
    limit = max_cov * THRESHOLD_PERCENT
    print(f"   Filtering nodes > {limit:.6f}...")
    indices = np.where(cov_col > limit)[0]
    n_sig = len(indices)
    print(f"   Nodes in bubble: {n_sig}")
    if n_sig == 0: sys.exit("Error: No significant covariance.")

    # 4. RECONSTRUCT FACTORS
    print("Reconstructing components...")
    
    # Distance
    dist_file = os.path.join(WDIST_DIR, f"station{real_station_id}.nc")
    if os.path.exists(dist_file):
        with Dataset(dist_file, 'r') as nc:
            d_all = nc.variables['dis'][:]
            dist_factor = np.exp(-(d_all[indices] / 111.0) / Lr)
    else:
        dist_factor = np.ones(n_sig)

    # Ratio
    if os.path.exists(DM_FILE):
        with Dataset(DM_FILE, 'r') as nc:
            dm_all = nc.variables['z'][:]
            dm_stn = dm_all[stn_node_idx] if dm_all[stn_node_idx] != 0 else 0.0001
            ratio = np.abs(dm_all[indices] / dm_stn)
            ratio_factor = np.where(ratio > 1, 1.0/ratio, ratio)
    else:
        ratio_factor = np.ones(n_sig)

    # Correlation
    r_factor = np.ones(n_sig)
    if os.path.exists(KCOEF_FILE):
        with Dataset(KCOEF_FILE, 'r') as nc:
            if 'station' in nc.variables:
                matches = np.where(nc.variables['station'][:].astype(int) == real_station_id)[0]
                if len(matches) > 0:
                    r_factor = nc.variables['coef'][indices, matches[0]]

    # 5. GENERATE MAPS
    df_mesh = read_fort14_coords(FORT14_FILE)
    df_bubble = df_mesh.iloc[indices].copy()
    
    df_bubble['cov'] = cov_col[indices]
    df_bubble['dist_fac'] = dist_factor
    df_bubble['ratio_fac'] = ratio_factor
    df_bubble['r_fac'] = r_factor
    
    lat = df_mesh.iloc[stn_node_idx]['lat']
    lon = df_mesh.iloc[stn_node_idx]['lon']

    print("\n--- Saving Maps (Jet + White SVG Text) ---")
    generate_map(df_bubble, 'cov', lat, lon, "map_1_COVAR.html", f"Total Cov (m_err^2={m_error_sq:.4f})")
    generate_map(df_bubble, 'dist_fac', lat, lon, "map_2_DIST.html", "Distance Factor")
    generate_map(df_bubble, 'ratio_fac', lat, lon, "map_3_RATIO.html", "Ratio Similarity")
    generate_map(df_bubble, 'r_fac', lat, lon, "map_4_CORR.html", "Correlation (r)")

    print("\n✅ Done. Open maps to verify white text.")

if __name__ == "__main__":
    main()
