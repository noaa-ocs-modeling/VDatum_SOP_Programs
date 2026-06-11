#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python3
import os
from pathlib import Path

def generate_vpop_in_files():
    # 1. Base Configuration Paths
    path_out = Path('out_vpop29_in_files')
    path_bp = Path('/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/output_BP2026_v7')
    path_marine = Path('/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/vgrid/out_marine_grid_nc')
    
    # Combined NetCDF source folder containing files like 'd_dd_diaPA_mhhw_Pac_SAL.nc'
    nc_source_dir = "/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/out_nc"
    
    # Tidal datums to step through
    varlist = ['mhhw', 'mhw', 'mlw', 'mllw', 'mtl', 'dtl']
    
    # Modes: Directly controls the suffix attached to the output files inside the .in file
    modes = [
        {"suffix": "",     "label": "standard datum (d)"},
        {"suffix": "_svu", "label": "SVU correction (diaPA)"}
    ]

    # Ensure the output directory exists
    path_out.mkdir(parents=True, exist_ok=True)

    if not path_bp.exists():
        print(f"Error: Bounding Polygon directory '{path_bp}' not found.")
        return
        
    # Gather all active bounding polygon directories
    bp_dirs = sorted([d for d in path_bp.iterdir() if d.is_dir()])
    print(f"Found {len(bp_dirs)} regional boundary directories to process.")

    for bp_path in bp_dirs:
        iname = bp_path.name
        print(f"--------------\nGenerating Templates For: {iname}")
        
        # Create separate .in files for standard datum and SVU variant passes
        for mode in modes:
            suffix = mode["suffix"]
            
            # 2. Build the configuration block string structure (NetCDF-to-NetCDF)
            in_content = f"""Input file for vpop29.f: Pacific Domain - {iname} ({mode['label']})
0  *** SHOULD BE ZERO *** <icorr_ord  (correct for tidal order)
0                           <icorr_dat  (correct to tidal datums)
0                           <ivoid_stop (stop if any pts unfilled
0                           <ismooth (number of times)
0  72 292  46 65     1389 86        <nprt, i,jprnt
0  0  0  0  0               <print:isv,itide,dv(# digits),ipond,mfill
-0.1524  0.0   lwd.nc       <value1,value2,new datum file name
none                        <MSF file
none                        <non-tidal
../../{path_bp.name}/{iname}/dpolygon_xyij01.dat                   <BP
none                                                    <grid bound polygon 
none                                                    <MSL NetCDF file
{iname}/unfilled                                         <output:unfilled grids
{iname}/itide                                            <output:tide points
{path_marine}/marine_{iname}.nc                        <marine grid
{len(varlist)}  1                        <nfmax(number of sets), ngrids
"""
            # 3. Inject NetCDF input paths and suffix-explicit NetCDF output paths
            for var in varlist:
                input_nc_file = f"{nc_source_dir}/d_dd_diaPA_{var}_Pac_SAL.nc"
                
                # Standard mode writes to 'PA_NE_HI_01/mhhw.nc'
                # SVU mode writes directly to 'PA_NE_HI_01/mhhw_svu.nc'
                output_nc_path = f"{iname}/{var}{suffix}.nc"
                
                in_content += f"{input_nc_file}\n"
                in_content += f"{output_nc_path}                     <output vdatum file name\n"
            
            # Ending block boundary flags
            in_content += "0  0\n"
            in_content += "--------------------------------------------------------------------------- \n"

            # 4. Save purely as a clean configuration file
            out_file = path_out / f"vpop29_{iname}{suffix}.in"
            with open(out_file, 'w') as f:
                f.write(in_content)
            print(f"  Wrote: {out_file}")

if __name__ == "__main__":
    generate_vpop_in_files()


# In[ ]:




