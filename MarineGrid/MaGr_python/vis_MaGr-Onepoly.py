import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import binary_dilation
import pandas as pd

# 1. Load the structured grid
# Path to your NetCDF file (Update this for Hawaii or CNMI accordingly)
#nc_path = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/marine_PA_W_CNMI_01.nc"
nc_path = "/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/marine_PA_NE_HI_01.nc"
print(f"Reading {nc_path}...")
with Dataset(nc_path, "r") as ds:
    field = ds.variables["field"][:]
    # box contains [alat0, alon0, dely, delx]
    alat0, alon0, dely, delx = ds.variables["box"][:]

jmax, imax = field.shape

# Calculate raw coordinates
clon_raw = alon0 + delx * np.arange(imax, dtype=np.float64)
clat = alat0 + dely * np.arange(jmax, dtype=np.float64)

# --- COORDINATE NORMALIZATION FOR QGIS ---
# This step ensures that Hawaii (often ~205E) is converted to (~ -155W)
# CNMI (~145E) will remain unchanged.
clon = np.where(clon_raw > 180, clon_raw - 360, clon_raw)
# -----------------------------------------

print("Applying Coastal Fringe filter (removing deep offshore points)...")

# 2. THE FRINGE FILTER
# Identify the expansion layers (Value 2 and above represent the inland/buffer zones)
expansion_mask = (field >= 2)

# Identify ONLY the water (Value 1) that is near the expansion
# iterations=10 creates a buffer to keep the immediate shoreline visible in QGIS
near_shore_water_area = binary_dilation(expansion_mask, iterations=10)
water_fringe_mask = (field == 1) & near_shore_water_area

# Combine: Expansion Layers + Near-shore Water Fringe
final_mask = expansion_mask | water_fringe_mask

# Find indices where our mask is True
j_idx, i_idx = np.where(final_mask)

if len(j_idx) > 0:
    # 3. Create DataFrame
    print(f"Extracting {len(j_idx)} points...")
    df = pd.DataFrame({
        'longitude': clon[i_idx],
        'latitude': clat[j_idx],
        'grid_value': field[j_idx, i_idx]
    })

    # 4. Save to CSV
    # Using a generic name or dynamic name based on the input
    #output_file = "Marine_CNMI_V8.csv"
    output_file = "./marine_HI_V9.csv"
    df.to_csv(output_file, index=False)
    
    print(f"\nSUCCESS: Created {output_file}")
    print(f"Longitude range in file: {df['longitude'].min()} to {df['longitude'].max()}")
    print("Action: Drag this CSV into QGIS. Set Geometry CRS to EPSG:4326.")
else:
    print("Error: No points found. Verify your NetCDF file contains values 1-11.")