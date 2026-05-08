#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import netCDF4 as nc
import numpy as np
import glob
import os

# --- CONFIGURATION ---
# Path to the folder containing your marine files
search_path = "./marine_*.nc" 

# List all matching files
files = glob.glob(search_path)
files.sort() # Keep them in order

print(f"Found {len(files)} files. Starting Audit...\n")
print(f"{'Filename':<35} | {'Unique Values':<25} | {'Status'}")
print("-" * 80)

bad_files = []

for nc_path in files:
    try:
        ds = nc.Dataset(nc_path, "r")
        # Read the field variable
        field = ds.variables['field'][:]
        unique_vals = np.unique(field)
        ds.close()

        # Check if it failed the expansion (only contains 0 and 1)
        # Note: We check if max value is <= 1
        if np.max(unique_vals) <= 1:
            status = "❌ FAILED (Only 0/1)"
            bad_files.append(nc_path)
        else:
            status = "✅ OK (Expansion Present)"

        filename = os.path.basename(nc_path)
        print(f"{filename:<35} | {str(unique_vals):<25} | {status}")

    except Exception as e:
        print(f"Error reading {nc_path}: {e}")

# --- SUMMARY ---
print("-" * 80)
print(f"Audit Complete.")
print(f"Total Files Checked: {len(files)}")
print(f"Files Needing Fix:   {len(bad_files)}")

if bad_files:
    print("\nList of files to re-run vgridder for:")
    for f in bad_files:
        print(f" - {f}")

