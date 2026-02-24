#!/usr/bin/env python3
"""
svu_export_tad_only.py
1. Reads Mesh (fort.14) for geometry.
2. Reads ONLY the TAD.xx file for data.
3. Exports a QGIS-ready UGRID NetCDF.
"""

import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# ================= CONFIGURATION =================
RUNID = "Pac_SAL"
DATUM = "mllw"  # Change this to extract different datums (e.g., 'msl', 'mllw')

# PATHS
BASE_DIR = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/"
FORT14 = os.path.join(BASE_DIR, "indata/fort.14")
TAD_FILE = "/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_01262026_TAD/TADxx.nc"

# OUTPUT FILE
OUT_FILE = os.path.join(BASE_DIR, f"qgis_{DATUM}_TAD_only.nc")

# TAD Index Mapping (Based on your provided script)
TAD_INDEX_MAP = {
    'mhhw': 0, 'mhw': 1, 'msl': 2, 'mtl': 3, 
    'dtl': 4, 'mlw': 5, 'mllw': 6
}
# =================================================

def save_qgis_file(filename, lat, triangles, lon_fixed, data_dict):
    """
    Writes UGRID file using the EXACT structure required for QGIS mesh layers.
    """
    print(f"   Writing: {os.path.basename(filename)}...")
    NP = len(lon_fixed)
    NE = len(triangles)

    with Dataset(filename, 'w', format='NETCDF4') as nc:
        # Dimensions
        nc.createDimension('node', NP)
        nc.createDimension('face', NE)
        nc.createDimension('nmax_face', 3)

        # Mesh Topology
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

        # Data Variables
        for var_name, var_data in data_dict.items():
            print(f"      + Layer: {var_name}")
            v = nc.createVariable(var_name, 'f4', ('node',), fill_value=-99999)
            v.mesh = 'mesh2d' 
            v[:] = var_data

def main():
    print(f"--- EXPORT TAD ONLY FOR {DATUM.upper()} ---")

    # 1. READ MESH (Required for QGIS to build the triangles)
    if not os.path.exists(FORT14):
        print(f"Error: Mesh file not found at {FORT14}")
        return

    print(f"Reading mesh: {os.path.basename(FORT14)}")
    with open(FORT14, 'r') as f:
        line1 = f.readline()
        line2 = f.readline().split()
        NE, NP = int(line2[0]), int(line2[1])

    # Efficiently read nodes and elements using pandas
    df_nodes = pd.read_csv(FORT14, delim_whitespace=True, header=None, 
                           skiprows=2, nrows=NP, usecols=[0,1,2], names=['id','x','y'])
    df_elems = pd.read_csv(FORT14, delim_whitespace=True, header=None, 
                           skiprows=2+NP, nrows=NE, usecols=[2,3,4], names=['n1','n2','n3'])
    
    # Prepare Mesh Data
    triangles = df_elems[['n1', 'n2', 'n3']].values - 1
    lon_nodes = df_nodes['x'].values
    lat_nodes = df_nodes['y'].values
    
    # GIS Fix (-180/180)
    lon_fixed = np.where(lon_nodes > 180, lon_nodes - 360, lon_nodes)

    # 2. READ TAD DATA
    print(f"Reading TAD file: {os.path.basename(TAD_FILE)}")
    tad_idx = TAD_INDEX_MAP.get(DATUM.lower())
    
    if tad_idx is None:
        print(f"Error: Datum '{DATUM}' not found in index map.")
        return

    tad_data = np.zeros(NP) # Default to zeros if read fails

    if os.path.exists(TAD_FILE):
        with Dataset(TAD_FILE, 'r') as nc:
            try:
                # Read specific datum row based on index
                raw_tad = nc.variables['datums'][tad_idx, :]
                
                # Safety check for size mismatch
                if len(raw_tad) >= NP:
                    tad_data = raw_tad[:NP]
                else:
                    print(f"Warning: TAD data length ({len(raw_tad)}) < Mesh Nodes ({NP}). Padding with zeros.")
                    tad_data[:len(raw_tad)] = raw_tad
            except KeyError:
                print("Error: Variable 'datums' not found in TAD file.")
            except Exception as e:
                print(f"Error reading TAD: {e}")
    else:
        print(f"Error: TAD file not found at {TAD_FILE}")
        return

    # 3. WRITE OUTPUT
    out_dict = {
        f'{DATUM}_tad': tad_data
    }
    
    save_qgis_file(OUT_FILE, lat_nodes, triangles, lon_fixed, out_dict)

    print("✅ Done! TAD file generated.")

if __name__ == "__main__":
    main()
