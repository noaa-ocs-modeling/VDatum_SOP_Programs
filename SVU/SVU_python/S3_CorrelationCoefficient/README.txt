
# Tidal Datum Correlation Calculator (Multi-Station / Parallel)

## Overview

The python code aligns with MATLAB-based logic for grid generation and spline fitting but is optimized for high-performance computing (HPC) using Python's `multiprocessing` and SLURM job arrays.

## 1\. Prerequisites & Environment

You need a Python 3 environment with the following libraries installed:

  * `numpy`
  * `xarray`
  * `scipy`
  * `pandas`
  * `netCDF4`

**Activate your environment:**

```bash
source /path/to/conda/bin/activate <your_env_name>
# Example: source /home/mojganr/miniconda3/bin/activate TAD
```

## 2\. Input Files

Ensure you have the absolute paths to the following files:

1. `TADhh.nc`: NetCDF file containing High Water times and values (`Htime_Hval`).
2.  `TADll.nc`: NetCDF file containing Low Water times and values (`Ltime_Lval`).
3.  `gageinfo.csv`: CSV file containing station metadata.
     Format: Must contain `station_id` (col 0) and `node_id` (col 1).

## 3\. Running the Analysis (SLURM)

There are two ways to run the analysis: processing **all datums** at once or running **individual datums**.

### Option A: Run All Datums (Recommended)

Use the master submission script `submit_coef_array_mp_3_multistations_All.sh`. This script loops through all datums (mlw, mllw, mhw, mhhw, mtl, dtl) and submits job arrays for each.

Before submitting, edit `submit_coef_array_mp_3_multistations_All.sh`:

1.  Set Paths: Update `HH_NC`, `LL_NC`, and `GAGE_CSV` to your input file locations.
2.  Set Config: Update `NCHUNKS` (total shards), `WORKERS` (CPUs per task), and `MAX_STATIONS`.
3.  Submit:
    ```bash
    sbatch submit_coef_array_mp_3_multistations_All.sh
    ```

### Option B: Run Individual Datums

If you only need to process one specific datum (e.g., MHHW), use the specific `.sh` file (e.g., `submit_coef_array_mp_3_multistations_mhhw.sh`).

1.  Open the specific `.sh` file.
2.  Update the input paths (`HH_NC`, `GAGE_CSV`) and configuration variables.
3.  Submit:
    ```bash
    sbatch submit_coef_array_mp_3_multistations_mhhw.sh
    ```

-----

## 4\. Post-Processing: Combining Output Files

The parallel jobs produce multiple "shard" files (e.g., `..._SLICE_0.nc`, `..._SLICE_1.nc`). You must combine these into a single NetCDF file using `combine_chucks.py`.

### IMPORTANT: Instructions for `combine_chucks.py`

The combination script is currently hardcoded for **MHHW**. To combine files for other datums (MLW, MHW, etc.), you must edit 3 specific lines in `combine_chucks.py` before running it.

#### Lines to Change:

  * Line 19 (Directory): Change the folder name where your shard files are located.
  * Line 20 (Input Pattern): Change the filename pattern to match the datum you are processing.
  * Line 47 (Output Filename): Change the name of the final merged file.

#### Example: Changing from MHHW to MLW

Open `combine_chucks.py` and modify:

Line 19:

  Change: `outdir = "out_mhhw_matlabAligned_multi"`
  To: `outdir = "out_mlw_matlabAligned_multi"`

Line 20:

  Change:`pattern = "coef_full_mhhw_matlabAligned_multi_SLICE_*.nc"`
  To: `pattern = "coef_full_mlw_matlabAligned_multi_SLICE_*.nc"`

Line 47:

  Change: `out_nc = "coef_full_mhhw_matlabAligned_multi.nc"`
  To: `out_nc = "coef_full_mlw_matlabAligned_multi.nc"`

#### How to Run:

After making the edits for your specific datum:

```bash
python3 combine_chucks.py
```

-----

## 5\. Output Format

The final combined NetCDF file will contain:

  Dimensions: `node` (total nodes processed), `station` (number of stations).
  Variable: `coef(node, station)` - The correlation coefficient (float32).
  Coordinates:
      * `node`: The absolute node IDs.
      * `station`: The station IDs from the CSV.


