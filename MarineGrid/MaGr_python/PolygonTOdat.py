#!/usr/bin/env python
# coding: utf-8


import os
from pathlib import Path
import shapefile

# ============================================================
# SETTINGS
# ============================================================
# Using the Meter-based version to match your resolution goals
INPUT_SHP = "/mnt/c/Users/mojgan.rostaminia/Documents/Hawaii_Pacific/MarinGrid/MarinGrid_polygons_4326_v7.shp"
OUTPUT_ROOT = Path("output_BP2026_v7")
WRAP_TO_360 = True  

# ============================================================
# HELPERS
# ============================================================
def wrap_lon(x):
    if not WRAP_TO_360: return x
    return x if x >= 0 else x + 360.0

def write_dat_files(folder_path, shp_geom):
    pts = shp_geom.points
    parts = list(shp_geom.parts) + [len(pts)]
    
    # We only care about the Primary/Exclusion Flag logic here
    # c: used for the main grid resolution boundaries
    configs = [
        ('c', True, False),  
        ('b', False, False), 
        ('d', False, True)   
    ]

    for suffix, use_flag, is_swapped in configs:
        filename = f"{suffix}polygon_xyij01.dat"
        with open(folder_path / filename, 'w') as f:
            
            for ring_idx in range(len(parts)-1):
                start_i = parts[ring_idx]
                end_i = parts[ring_idx+1]
                ring_coords = pts[start_i:end_i]
                
                # REVISED: No Flag 2 allowed anywhere.
                # Every ring is treated as a Primary Boundary (Flag 0).
                # Step 3 Shoreline will handle the exclusion of land cells.
                outer_flag, ring_id = 0, 1
                
                col4 = outer_flag if use_flag else ring_id

                for i, (x, y) in enumerate(ring_coords):
                    xw = wrap_lon(x)
                    # 1 for ring closure, 0 for internal nodes
                    is_closed = 1 if (i == 0 or i == len(ring_coords) - 1) else 0
                    
                    v1 = y if is_swapped else xw
                    v2 = xw if is_swapped else y
                    
                    # Consistent 12.6f formatting for coordinate precision
                    f.write(f"{v1:12.6f} {v2:12.6f} {is_closed:3d} {col4:6d}\n")

# ============================================================
# MAIN
# ============================================================
def main():
    if not OUTPUT_ROOT.exists():
        OUTPUT_ROOT.mkdir(parents=True)

    print(f"Reading Shapefile: {INPUT_SHP}")
    sf = shapefile.Reader(INPUT_SHP)
    records = sf.shapeRecords()

    for sr in records:
        # Clean naming for folder structures
        name = sr.record["MATNAME"].strip().replace('-', '_')
        folder_path = OUTPUT_ROOT / f"{name}_01"
        folder_path.mkdir(exist_ok=True)

        print(f"Exporting: {name}_01 with all Flag 0")
        write_dat_files(folder_path, sr.shape)

if __name__ == "__main__":
    main()





