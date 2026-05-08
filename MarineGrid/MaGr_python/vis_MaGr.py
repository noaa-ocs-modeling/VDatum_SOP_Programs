import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import binary_dilation
import pandas as pd
import glob
import os

# --- 1. SETUP PATHS ---
input_dir = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/"
output_dir = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/outputs/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# FIXED LINE: This correctly searches for files starting with 'marine'
nc_files = glob.glob(os.path.join(input_dir, "marine*.nc"))

print(f"Starting process for {len(nc_files)} files...")

# --- 2. LOOP OVER ALL FILES ---
for nc_path in nc_files:
    base_filename = os.path.basename(nc_path)
    output_file = os.path.join(output_dir, base_filename.replace(".nc", ".csv"))
    
    print(f"Reading {nc_path}...")
    
    try:
        with Dataset(nc_path, "r") as ds:
            field = ds.variables["field"][:]
            alat0, alon0, dely, delx = ds.variables["box"][:]

        jmax, imax = field.shape

        clon_raw = alon0 + delx * np.arange(imax, dtype=np.float64)
        clat = alat0 + dely * np.arange(jmax, dtype=np.float64)

        # --- COORDINATE NORMALIZATION FOR QGIS ---
        clon = np.where(clon_raw > 180, clon_raw - 360, clon_raw)

        # --- 3. THE FRINGE FILTER ---
        expansion_mask = (field >= 2)
        near_shore_water_area = binary_dilation(expansion_mask, iterations=10)
        water_fringe_mask = (field == 1) & near_shore_water_area
        final_mask = expansion_mask | water_fringe_mask

        j_idx, i_idx = np.where(final_mask)

        if len(j_idx) > 0:
            # 4. Create DataFrame
            df = pd.DataFrame({
                'longitude': clon[i_idx],
                'latitude': clat[j_idx],
                'grid_value': field[j_idx, i_idx]
            })

            # 5. Save to CSV
            df.to_csv(output_file, index=False)
            print(f"SUCCESS: Created {output_file} ({len(j_idx)} points)")
        else:
            print(f"No points found in {base_filename}. Skipping.")

    except Exception as e:
        print(f"FAILED to process {base_filename}: {e}")

print("\nFinished! All matching 'marine' files have been processed.")

print("\nFinished! All files have been processed individually.")