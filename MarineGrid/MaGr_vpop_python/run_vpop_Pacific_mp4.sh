#!/bin/bash --login
#SBATCH -J vpop_Pacific
#SBATCH -t 6:00:00
#SBATCH --array=0-48
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -o logs/vpop29_PA_%A_%a.out
#SBATCH -e logs/vpop29_PA_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov

set -euo pipefail
set -x

# ------------------------------------------------------------
# Modules / Python environment
# ------------------------------------------------------------
source /home/mojganr/miniconda3/bin/activate TAD

# Avoid oversubscription (same as vgridder)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ------------------------------------------------------------
# User settings - paths for vpop
# ------------------------------------------------------------
# Directory containing vpop29 input files (.in)
export VPOP_INDIR="/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/vpop/out_vpop29_in_files"

# Directory containing marine grid NetCDF files (from vgridder)
export MARINE_GRID_DIR="/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/vgrid/out_marine_grid_nc"

# Output directory
export VPOP_OUTDIR="/work2/noaa/vdatum/mojganr/work_adcirc/MaGr/vpop/out_vpop29_files_nc"

# Shared Pacific ADCIRC mesh needed because datum .nc files contain values
# but not mesh geometry/connectivity
export MESH_PATH="/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/PA_6sec_SAL_08182025/fort.14"

#set mode explicitly
#export VPOP_MODE="standard"
export VPOP_MODE="svu"

# Number of polygons (0-48 = 49 total)
export NPOLY=49

# Workers for any internal parallelism
export VPOP_NWORKERS="${SLURM_CPUS_PER_TASK:-4}"

# ------------------------------------------------------------
# Create logs directory and run
# ------------------------------------------------------------
mkdir -p logs
mkdir -p ${VPOP_OUTDIR}

python -u vpop_nc_Pacific_mp4.py

echo "Done with task ${SLURM_ARRAY_TASK_ID}"
