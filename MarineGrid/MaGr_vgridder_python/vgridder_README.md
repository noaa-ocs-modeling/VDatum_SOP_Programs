# VDatum Marine Grid Generation Pipeline (`vgridder`)

## Overview
This repository contains the complete pipeline for generating high-resolution, structured marine grids for the VDatum project. The workflow transitions from local GIS preprocessing (SMS/QGIS) to automated, massively parallel grid generation on the NOAA Hercules HPC cluster. 

The core engine (vgridder_nc_Pacific_mp4.py) utilizes vectorized boundary enforcement and bounding-box pre-filtering to rapidly process tens of millions of grid points while properly handling nested island domains (holes). Edge cases, such as riverine discontinuities or buffer zone protections, are managed strictly through a localized command-file architecture.

The legacy Fortran vgridder architecture relied on brute-force, point-by-point nested loops to calculate boundary intersections and layer expansions, resulting in severe computational bottlenecks for massive domains like the 500m Pacific Ocean grid. During the modernization to Python, the core logic was fundamentally optimized rather than just translated. By replacing the iterative pixel-by-pixel checks with NumPy vectorization (np.searchsorted) and introducing geographic bounding-box pre-filters, the new Python engine evaluates entire latitude rows simultaneously while instantly skipping irrelevant masked regions. This architectural shift leverages highly optimized C-backed array operations to deliver massive performance gains. Crucially, this speed does not compromise accuracy; when benchmarked using the complex Alaska regional domain, the new Python code reduced the processing time from several hours down to just 30 minutes, while producing grids that were 100% identical to the legacy Fortran outputs, proving strict mathematical equivalence.

***
---

## Script Inventory

Phase	Script Name	                        Location	 Primary Function
QA/QC	Test_ovlappoly_Pacific.py	        Desktop	     Validates internal polygon alignment and prevents overlaps.
QA/QC	Test_ovlappoly_Pacific&Alaska.py	Desktop	     Validates external boundary alignment against neighboring project grids.
Prep	PolygonTOdat.py	                    Desktop	     Converts Shapefiles to the required VDatum .dat boundary format.
Config	vgridder_inGenerator.py	            Hercules     Generates the .in configuration files for the gridder engine.
Compute	run_vgridder_Pacific_mp4.sh	        Hercules   	 SLURM batch execution script for the gridder array.
Compute	vgridder_nc_Pacific_mp4.py	        Hercules     The core Python grid generation engine.
Audit	Check_marine_nc-onepoly.py	        Hercules     Prints the unique array values of a single NetCDF to verify the exact number of buffer layers generated (e.g., 0, 1, 2... 10).
Audit	Check_marine_nc.py	                Hercules     Audits all generated NetCDF files to ensure layer generation succeeded globally.
Visual	vis_MaGr-Onepoly.py	                Hercules 	 Converts a single NetCDF to a sub-sampled CSV for QGIS visualization.
Visual	vis_MaGr.py	                        Hercules     Converts all NetCDFs to CSVs for global QGIS visualization.
---

## Phase 1: Local Pre-Processing (Desktop)
Before moving to the HPC, the raw bounding polygons must be generated, validated, and converted into the correct text formats.

1. **Draw and Clean Polygons:** Generate initial bounding polygons in SMS, then import into QGIS to clean geometries, name features, and verify boundaries.
2. **Internal Alignment Check:** Run `Test_ovlappoly_Pacific.py` to ensure none of the internal project polygons overlap or cross into one another. **Required CRS: EPSG 3832.**
3. **External and Internal Alignment Check:** Run `Test_ovlappoly_Pacific&Alaska.py` to guarantee the outer project boundary seamlessly aligns with neighboring VDatum regions (e.g., Alaska, West Coast). **Required CRS: EPSG 3832.**
4. **Generate Boundary Files:** Run `PolygonTOdat.py` to convert the validated shapefile into individual VDatum boundary files (`cpolygon_xyij01.dat`). **Required CRS: EPSG 4326.** 
5. **Transfer:** Upload the resulting `output_BP` folder to the Hercules cluster.

---

## Phase 2: Configuration & Execution (Hercules HPC)

1. **Generate Configurations:** Run `vgridder_inGenerator.py`. This script dynamically builds the `vgridder_[name].in` files, automatically applying high-resolution settings (e.g., 0.00045 degrees) to island domains and low-resolution settings (e.g., 0.002 degrees) to the massive `PA_C_Ocean_01` domain.
2. **Execute SLURM Array:** Submit the job to the cluster using `sbatch run_vgridder_Pacific_mp4.sh`. This triggers `vgridder_nc_Pacific_mp4.py` across all polygons simultaneously using parallel processing.

### Edge Case Management (`local_overrides.dat`)
For specific geographic anomalies (e.g., severing artificial land bridges or forcing disconnected rivers to stay wet), do not alter the global Python script. Instead, place a `local_overrides.dat` file inside the specific polygon's boundary folder (e.g., `PA_NE_HI_01/local_overrides.dat`).
* **Format:** `COMMAND LON_MIN LON_MAX LAT_MIN LAT_MAX`
* **Available Commands:**
    * `FORCE_LAND`: Converts the specified bounding box to land (`0`), forcing the layer engine to treat it as a hard barrier.
    * `FORCE_WATER`: Converts the specified bounding box to water (`1`) and actively protects it from the inland pond-eraser logic.
* **Example:** `FORCE_WATER -155.85 -155.84 20.01 20.02`

---

## Phase 3: QA/QC and Visualization
The resulting NetCDF grids must be audited for successful buffer layer expansion and visually verified before proceeding to the `vpop` stage.

1. **Single File Layer Audit:** Run `Check_marine_nc-onepoly.py` to print the unique values of a specific output grid. This confirms exactly which inland buffer layers (2, 3, 4, etc.) were successfully generated.
2. **Global Audit:** Run `Check_marine_nc.py` to scan all generated marine NetCDF files. This will flag any grid that failed to generate buffer layers (i.e., files that only contain `0` and `1`).
3. **Single File Visualization:** Run `vis_MaGr-Onepoly.py` to extract a single marin.nc file for each polygon into a CSV format. This applies a Fringe Filter (or Geographic Bounding Box) and a Sub-sampling Stride to keep the file size safe for desktop GIS memory limits.
4. **Global Visualization:** Run `vis_MaGr.py` to bulk-process all successfully generated NetCDFs into visualizable CSVs. 
5. **Final Check:** Drag the resulting CSVs into QGIS (Set CRS to EPSG:4326) and apply Categorized Symbology to the `grid_value` field to visually verify the boundaries and buffer layers.
