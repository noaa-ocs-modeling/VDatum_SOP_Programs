#!/bin/bash
#SBATCH -J vgrid_Pac_v2mp4
#SBATCH -t 1:00:00
# --- CHANGE THIS TO TEST ONLY TWO ISLANDS ---
# Example: Use the indices for Hawaii and another island (e.g., 4 and 15)
#SBATCH --array=1,44
# ---------------------------------------------
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH -o logs/vgrid_Pac_v2mp4_%A_%a.out
#SBATCH -e logs/vgrid_Pac_v2mp4_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov


set -euo pipefail
source /home/mojganr/miniconda3/bin/activate TAD

export NPOLY=49
export BP_ROOT="/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/output_BP2026_v7"
export COAST_PATH="/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/Pacific_shoreline_merged_nozm_wgs84.shp"
export MESH_PATH="/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/PA_6sec_SAL_08182025/fort.14"
export DATUM_FILE="/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/out_nc/d_dd_diaPA_mllw_Pac_SAL.nc"
export REF_DIR="/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/vgrid/out_marine_grid_in_files/"

mkdir -p logs
export VGRID_NWORKERS=${SLURM_CPUS_PER_TASK:-4}

python -u vgridder_v24nc4_Pacific_v12_mp4.py
