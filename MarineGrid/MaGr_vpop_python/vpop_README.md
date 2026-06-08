## Overview

This repository contains the complete pipeline for the VPOP phase of the VDatum project. This workflow transitions the legacy tidal datum interpolation process from a sequential Fortran architecture into an automated, massively parallelized Python pipeline executed on the NOAA Hercules HPC cluster.

The core engine (`vpop_nc_Pacific_mp4.py`) takes the structured marine grids and populates them with corrected value of tidal datums (e.g., MHHW, MLLW) and SVU (Spatially Varying Uncertainty).

The legacy Fortran `vpop29nc4.f` architecture relied on brute-force nearest-node searches and point-by-point iterative sweeps to perform triangle interpolations and neighbor-weighted fills. During the modernization to Python, the engine was fundamentally accelerated without losing mathematical parity. By integrating `scipy.spatial.cKDTree` for spatial indexing, the script instantly identifies candidate ADCIRC nodes for interpolation rather than scanning the entire mesh. Furthermore, critical bottlenecks like the Fortran-synchronous neighbor-weighted fill (Step 7) and planar triangle interpolation (Step 5) were reconstructed using `Numba` JIT (Just-In-Time) compilation. This guarantees 100% mathematical equivalence to the original Fortran outputs while leveraging modern vectorization to process massive Pacific and Alaskan domains in a fraction of the time.

---

## Script Inventory

| Phase | Script Name | Location | Primary Function |
| --- | --- | --- | --- |
| Config | `vpop_inGenerator.py` | Hercules HPC | Generates the `.in` configuration files for standard and SVU correction runs. |
| Compute | `run_vgrid_Pacific_mp4.sh` | Hercules HPC | SLURM batch execution script (Note: Provided file acts as a template for array submission, analogous to the `run_vpop29` equivalent). |
| Compute | `vpop_nc_Pacific_mp4.py` | Hercules HPC | The core Python VPOP generation engine that populates marine grids with ADCIRC datums. |
| Visual | `export_vpop_structured_points_qgis_HI.py` | Hercules HPC / Local | Converts the populated NetCDF fields into QGIS-ready CSVs, applying a coastal fringe mask to manage file sizes. |

---

## Phase 1: Execution (Hercules HPC)
* **Generate Configurations:** Run vpop_inGenerator.py. This script dynamically builds the vpop29_[name].in files, stepping through standard tidal datums (MHHW, MHW, MLW, MLLW, MTL, DTL) and establishing output paths. It simultaneously creates secondary configuration blocks specifically for SVU (diaPA) uncertainty passes.
* **Execute SLURM Array:** Submit the job to the cluster using the appropriate SLURM array script (patterned after `run_vgrid_Pacific_mp4.sh`). This triggers the core Python engine across all regional polygons simultaneously using parallel processing, mapping each SLURM task ID to a specific `.in` file.

---

## Phase 2: The Core Vectorized Interpolation Engine

The `vpop_nc_Pacific_mp4.py` script replaces the legacy Fortran codebase. It executes a rigorous, multi-step fill process to guarantee data propagation from the unstructured ADCIRC mesh to the high-resolution structured marine grid.

* **Triangle Interpolation (Step 5):** Uses a `cKDTree` to identify the two closest usable ADCIRC nodes to any target grid cell. A Numba-accelerated Fortran port then performs exact planar interpolation (`inside_tri3` and `abc_eval`) to fill primary water cells.
* **Small-Scale Fill (Step 6):** Addresses gaps by applying a highly localized node-in-cell bounding box search with a configurable tolerance (`VPOP_SMALL_EPS`), taking a weighted average of nearby nodes.
* **Weighted Neighbor Fill (Step 7):** Iteratively fills remaining inland tidal cells. To ensure strict parity with Fortran, this step utilizes a Numba-compiled synchronous sweep (`_fill_weighted_fortran_sync_core`) that calculates temporary values across the grid and only updates the actual matrix after the sweep finishes.
* **Layer Fill (Step 8):** Pushes the tidal datum values up riverine systems and into buffer zones (layers 2+) generated during the `vgridder` phase.

---

## Phase 3: QA/QC and Visualization

Because standard 500m Pacific grids contain tens of millions of nodes, direct visualization of the entire VPOP field can crash local GIS software.

* **Coastal Fringe Extraction:** Run `export_vpop_structured_points_qgis_HI.py`. Instead of dumping the entire NetCDF, this script utilizes `scipy.ndimage.binary_dilation` to isolate only the structured grid points that reside along the critical land-water boundary (the coastal fringe).
* **Final Check:** Drag the resulting CSV into QGIS (Set Geometry CRS to EPSG:4326). Apply Graduated Colors to the `vpop_value` field to visually verify that tidal interpolations smoothly transition across the coastal fringe without anomalies or disconnected void zones.
