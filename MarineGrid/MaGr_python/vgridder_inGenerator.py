#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import math
import numpy as np

def generate_in_files():
    # 1. Configuration
    path_out = 'out_marine_grid_in_files'
    pathpc_bps = '/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/output_BP2026_v7'
    
    # Ensure the output directory exists
    os.makedirs(path_out, exist_ok=True)

    # Find all BP directories inside pathpc_bps
    if not os.path.exists(pathpc_bps):
        print(f"Error: Directory '{pathpc_bps}' not found.")
        return
        
    bp_dirs = [d for d in os.listdir(pathpc_bps) if os.path.isdir(os.path.join(pathpc_bps, d))]
    bp_dirs.sort()

    for iname in bp_dirs:
        print(f"--------------\nProcessing: {iname}")
        
        
        # 2. Pacific Domain Logic (Resolution and Inland Distance)
        ielim_ponds = 1
        # We define what is NOT high-res
        low_res_folders = ["PA_C_Ocean_01"]
        
        if not any(low_res in iname for low_res in low_res_folders):
            # Settings for high-res islands (50 meters)
            dx = 0.00045
            dy = 0.00045
            # 500m inland extension (Only applies to land/islands)
            dinland = 500.0 / 111120.0  
            nlayer = math.ceil(dinland / dy)
            print(f"  dx={dx:.4f}, dy={dy:.6f}, Inland: 500m -> Layers: {nlayer}")
        else:
            # Settings for the main/broad ocean (1 km)
            dx = 0.009
            dy = 0.009
            # The Ocean doesn't touch land, so it gets NO inland extension
            dinland = 0.0 
            nlayer = 0
            print(f"  dx={dx:.4f}, dy={dy:.6f}, Inland: 0m -> Layers: {nlayer}")        
        # 3. Read Bounding Polygon to calculate xy_lim
        poly_file = os.path.join(pathpc_bps, iname, 'cpolygon_xyij01.dat')
        if not os.path.exists(poly_file):
            print(f"  Warning: {poly_file} not found. Skipping.")
            continue
            
        # Read the DAT file to find the min/max coordinates
        data = np.genfromtxt(poly_file)
        min_x, max_x = np.min(data[:, 0]), np.max(data[:, 0])
        min_y, max_y = np.min(data[:, 1]), np.max(data[:, 1])
        
        # Calculate bounding box limits with padding
        alat1 = min_y - dy
        alat2 = max_y + dy
        alon1 = min_x - dx
        alon2 = max_x + dx
        
        # 4. Generate the exact .in file structure
        in_content = f"""file for vgridder_v24nc4_Pacific_v2_mp4.py: Pacific VDatum {pathpc_bps}/
4   {alat1:.6f}  {alat2:.6f} {alon1:.6f} {alon2:.6f}   < {iname}
{dx:.5f} {dy:.6f}                 <delx,dely (degrees)
{ielim_ponds}                           <ielim_ponds (of size <=ielim_ponds)
0                            <nadjust, i,j, jfield
{nlayer}                            <layers
0   0      (was 1  0)        <nobarr1,nobarr2
0  70 169                    <npr,iprt(),jprt()
1  1                         <numnulls=nulls to be considered dry. nblock 
/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/out_nc/d_dd_diaPA_mllw_Pac_SAL.nc   <1.ADCIRC datum, or none 
none                         <2.EEZ
/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/Pacific_shoreline_merged_nozm_wgs84.shp   <3.coastline
../{pathpc_bps}/{iname}/cpolygon_xyij01.dat                  <4.BP
marine_{iname}                     <5.output Marine Grid file in vpop.f format
marine_{iname}                     <6.output Marine Grid file in GTX format
water_{iname}                      <7.file with original water points (w/o main)
ponds_{iname}                      <8.file with ponds (w/o main)
dry_{iname}                        <9.file with land/water plus dried nodes (-1)
--------------------------------------------------------------------------- 
"""

        # 5. Write to File
        out_file = os.path.join(path_out, f"vgridder_{iname}.in")
        with open(out_file, 'w') as f:
            f.write(in_content)
        print(f"  Wrote: {out_file}")

if __name__ == "__main__":
    generate_in_files()


# Resolution,Layers,Visual Window,Physical Search
# 500 m, 0.0020° (Broad),3,7×7 cells,~600m
# 50 m, 0.0005° (Hawaii),3,7×7 cells,~150m
