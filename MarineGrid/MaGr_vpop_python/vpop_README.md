## Overview
This repository contains the complete pipeline for the VPOP phase of the VDatum project. The workflow modernizes the legacy tidal datum interpolation process from a sequential Fortran architecture into an automated, massively parallelized Python pipeline executed on the NOAA Hercules HPC cluster. The core engine, `vpop_nc_Pacific_mp4.py`, reads the structured marine grids and populates them with corrected tidal datum values (e.g., MHHW, MLLW) and SVU (Spatially Varying Uncertainty). The legacy Fortran `vpop29nc4.f` architecture relied on brute-force nearest-node searches and point-by-point iterative sweeps to perform triangle interpolation and neighbor-weighted filling. During the modernization to Python, the engine was fundamentally accelerated while preserving practical parity with the legacy workflow. By integrating `scipy.spatial.cKDTree` for spatial indexing, the script efficiently identifies candidate ADCIRC nodes for interpolation instead of scanning the entire mesh. In addition, critical bottlenecks such as the Fortran-synchronous neighbor-weighted fill (Step 7) and planar triangle interpolation (Step 5) were reconstructed using `Numba` JIT (Just-In-Time) compilation. This allows the workflow to retain the original computational logic while leveraging modern vectorization to process large multi-polygon coastal domains much more efficiently.

---

## Script Inventory
| Phase | Script Name | Location | Primary Function |
| --- | --- | --- | --- |
| Config | `vpop_inGenerator.py` | Hercules HPC | Generates the `.in` configuration files for standard and SVU correction runs. |
| Compute | `run_vpop_Pacific_mp4.sh` | Hercules HPC | SLURM batch execution script for array-based submission of polygon jobs. |
| Compute | `vpop_nc_Pacific_mp4.py` | Hercules HPC | Core Python VPOP engine that populates marine grids with interpolated ADCIRC datums. |
| Visual | `export_vpop_structured_points_qgis_HI.py` | Hercules HPC / Local | Converts populated NetCDF fields into QGIS-ready CSVs, applying a coastal fringe mask to reduce file size. |

---

## Phase 1: Execution (Hercules HPC)
* **Generate Configurations:** Run `vpop_inGenerator.py`. This script dynamically builds the `vpop29_[name].in` files, stepping through standard tidal datums (MHHW, MHW, MLW, MLLW, MTL, DTL) and establishing output paths. It also creates secondary configuration blocks for SVU (`diaPA`) uncertainty runs.
* **Execute SLURM Array:** Submit the job to the cluster using the appropriate SLURM array script, such as `run_vpop_Pacific_mp4.sh`. This launches the core Python engine across all regional polygons in parallel, mapping each `SLURM_ARRAY_TASK_ID` to a specific `.in` file.

---

## Phase 2: The Core Vectorized Interpolation Engine
The `vpop_nc_Pacific_mp4.py` script reimplements the legacy Fortran workflow in Python and executes a rigorous multi-step fill process to propagate values from the unstructured ADCIRC mesh to the high-resolution structured marine grid.

* **Triangle Interpolation (Step 5):** Uses a `cKDTree` to identify the two closest usable ADCIRC nodes to each target grid cell. A Numba-accelerated Fortran-style implementation then performs planar interpolation (`inside_tri3` and `abc_eval`) to fill primary water cells.

* **Small-Scale Fill (Step 6):** Addresses remaining gaps using a highly localized node-in-cell bounding-box search with a configurable tolerance (`VPOP_SMALL_EPS`), then applies a weighted average of nearby nodes.

* **Weighted Neighbor Fill (Step 7):** Iteratively fills remaining tidal cells. To preserve Fortran behavior, this step uses a Numba-compiled synchronous sweep (`_fill_weighted_fortran_sync_core`) that computes temporary values across the grid and updates the main arrays only after each sweep is complete.

* **Layer Fill (Step 8):** Extends tidal datum values into riverine systems and buffer zones (layers 2+) generated during the `vgridder` phase.
---

## Phase 3: QA/QC and Visualization
Because standard 500 m Pacific grids contain tens of millions of nodes, direct visualization of the full VPOP field can overwhelm local GIS software.

* **Coastal Fringe Extraction:** Run `export_vpop_structured_points_qgis_HI.py`. Instead of exporting the entire NetCDF, this script uses `scipy.ndimage.binary_dilation` to isolate only the structured grid points that lie along the critical land-water boundary (the coastal fringe).

* **Final Check:** Load the resulting CSV into QGIS and set the geometry CRS to `EPSG:4326`. Apply graduated colors to the `vpop_value` field to visually verify that tidal interpolations transition smoothly across the coastal fringe without obvious artifacts.hout anomalies or disconnected void zones.
