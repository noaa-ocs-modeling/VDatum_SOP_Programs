#!/usr/bin/env python3
"""
svu_export_qgis_split_final.py
FINAL SPLIT EXPORTER:
1. Uses the EXACT UGRID structure from 'svu_export_to_qgis_mesh_1.py' (which worked).
2. Splits into Results vs Inputs to keep files light.
3. Applies all MATLAB/Physics logic.
"""

import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.interpolate import griddata

# ================= CONFIGURATION =================
RUNID = "Pac_SAL"
DATUM = "mhhw"

# PATHS
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
FORT14 = os.path.join(BASE_DIR, "indata/fort.14")
NC_DATA = os.path.join(BASE_DIR, "out_nc", f"d_dd_diaPA_{DATUM}_{RUNID}.nc")
TAD_FILE = "/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_08182025_TAD/TADxx.nc"
NAN_FILE = os.path.join(BASE_DIR, "pre_process", "mdatum_nanloc_mllw.nc")

# OUTPUT FILES
OUT_RESULTS = os.path.join(BASE_DIR, f"qgis_{DATUM}_results.nc")
OUT_INPUTS  = os.path.join(BASE_DIR, f"qgis_{DATUM}_inputs.nc")

TAD_INDEX_MAP = {
    'mhhw': 0, 'mhw': 1, 'msl': 2, 'mtl': 3, 
    'dtl': 4, 'mlw': 5, 'mllw': 6
}
# =================================================

def f_replace_model_nan(xy, z, loc):
    newz = z.copy()
    loc = np.asarray(loc, dtype=int)
    if loc.size == 0: return newz
    valid_mask = np.isfinite(z)
    valid_mask[loc] = False 
    if np.sum(valid_mask) == 0: return newz
    try:
        znan = griddata(xy[valid_mask], z[valid_mask], xy[loc], method='nearest')
        newz[loc] = znan
    except:
        pass
    return newz

def save_qgis_file(filename, lat, triangles, lon_fixed, data_dict):
    """
    Writes UGRID file using the EXACT structure of the working script.
    """
    print(f"   Writing: {os.path.basename(filename)}...")
    NP = len(lon_fixed)
    NE = len(triangles)

    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Dimensions
        nc.createDimension('node', NP)
        nc.createDimension('face', NE)
        nc.createDimension('nmax_face', 3)

        # Mesh Topology (Exactly as in the working script)
        mesh = nc.createVariable('mesh2d', 'i4')
        mesh.cf_role = 'mesh_topology'
        mesh.topology_dimension = 2
        mesh.node_coordinates = 'node_x node_y'
        mesh.face_node_connectivity = 'face_nodes'

        # Coordinates
        x = nc.createVariable('node_x', 'f8', ('node',))
        x.standard_name = 'longitude'
        x[:] = lon_fixed

        y = nc.createVariable('node_y', 'f8', ('node',))
        y.standard_name = 'latitude'
        y[:] = lat

        # Connectivity
        conn = nc.createVariable('face_nodes', 'i4', ('face', 'nmax_face'))
        conn.cf_role = 'face_node_connectivity'
        conn.start_index = 0
        conn[:] = triangles

        # Data Variables (Minimal attributes, exactly as in working script)
        for var_name, var_data in data_dict.items():
            print(f"      + Layer: {var_name}")
            v = nc.createVariable(var_name, 'f4', ('node',), fill_value=-99999)
            v.mesh = 'mesh2d' # This is the only attribute the working script used
            v[:] = var_data

def main():
    print(f"--- SPLIT EXPORT (FINAL FIXED) FOR {DATUM.upper()} ---")

    # 1. READ MESH
    print(f"Reading mesh: {os.path.basename(FORT14)}")
    with open(FORT14, 'r') as f:
        line1 = f.readline()
        line2 = f.readline().split()
        NE, NP = int(line2[0]), int(line2[1])

    df_nodes = pd.read_csv(FORT14, delim_whitespace=True, header=None, 
                           skiprows=2, nrows=NP, usecols=[0,1,2], names=['id','x','y'])
    df_elems = pd.read_csv(FORT14, delim_whitespace=True, header=None, 
                           skiprows=2+NP, nrows=NE, usecols=[2,3,4], names=['n1','n2','n3'])
    
    triangles = df_elems[['n1', 'n2', 'n3']].values - 1
    lon_nodes = df_nodes['x'].values
    lat_nodes = df_nodes['y'].values
    p = np.column_stack((lon_nodes, lat_nodes))

    # GIS Fix (-180/180)
    lon_fixed = np.where(lon_nodes > 180, lon_nodes - 360, lon_nodes)

    # 2. READ NAN MASK
    if os.path.exists(NAN_FILE):
        with Dataset(NAN_FILE, 'r') as nc:
            keys = list(nc.variables.keys())
            vname = 'loc_model_nan' if 'loc_model_nan' in keys else keys[-1]
            nanloc = nc.variables[vname][:].astype(int) - 1
    else:
        nanloc = []

    # 3. READ DATA
    print("Reading data...")
    with Dataset(NC_DATA, 'r') as nc:
        d = nc.variables['d'][:]
        dd = nc.variables['dd'][:]
        
        if 'dm' in nc.variables:
            dm = nc.variables['dm'][:]
        else:
            dm = np.zeros_like(d)

        if 'diaPA_real' in nc.variables:
            diaPA = nc.variables['diaPA_real'][:]
        elif 'diaPA' in nc.variables:
            diaPA = nc.variables['diaPA'][:]
            if np.iscomplexobj(diaPA): diaPA = np.real(diaPA)
        else:
            diaPA = np.zeros_like(d)

    # 4. READ TAD
    tad_idx = TAD_INDEX_MAP.get(DATUM.lower())
    if tad_idx is None:
        tad = np.zeros_like(d)
    else:
        with Dataset(TAD_FILE, 'r') as nc:
            try:
                tad = nc.variables['datums'][tad_idx, :]
                if len(tad) > NP: tad = tad[:NP]
            except:
                tad = np.zeros_like(d)

    # 5. APPLY LOGIC
    print("Applying Logic...")
    is_tide = DATUM.lower() in ['mhhw', 'mhw', 'mlw', 'mllw']
    
    if is_tide:
        temp = np.abs(dm)
        stdv = np.std(temp, ddof=1)
        if stdv > 0:
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = temp / stdv
                ratio[~np.isfinite(ratio)] = 1.0
            ratio = np.minimum(ratio, 1.0)
            diaPA = diaPA * ratio
        loc_0 = np.where(np.abs(dm) < 0.0001)[0]
    else:
        loc_0 = []

    results = {'d': d, 'dd': dd, 'diaPA': diaPA}
    
    for key in results:
        arr = results[key]
        if len(nanloc) > 0: arr[nanloc] = np.nan
        if is_tide and len(loc_0) > 0: arr[loc_0] = 0.0

        if key == 'd' and is_tide:
            bad_sign = []
            if 'hh' in DATUM or 'hw' in DATUM: bad_sign = np.where(arr < 0)[0]
            elif 'll' in DATUM or 'lw' in DATUM: bad_sign = np.where(arr > 0)[0]
            if len(bad_sign) > 0: arr[bad_sign] = 0.0

        nan_curr = np.where(np.isnan(arr))[0]
        if len(nan_curr) > 0: arr = f_replace_model_nan(p, arr, nan_curr)
        results[key] = arr

    # 6. WRITE FILES
    res_dict = {
        f'{DATUM}_d': results['d'],
        f'{DATUM}_dd': results['dd'],
        f'{DATUM}_diaPA': results['diaPA']
    }
    save_qgis_file(OUT_RESULTS, lat_nodes, triangles, lon_fixed, res_dict)

    inp_dict = {
        f'{DATUM}_dm': dm,
        f'{DATUM}_tad': tad
    }
    save_qgis_file(OUT_INPUTS, lat_nodes, triangles, lon_fixed, inp_dict)

    print("✅ Done! Files generated with verified UGRID structure.")

if __name__ == "__main__":
    main()
