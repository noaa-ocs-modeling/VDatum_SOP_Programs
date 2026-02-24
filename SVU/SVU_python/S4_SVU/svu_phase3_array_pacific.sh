#!/bin/bash
#SBATCH -J svu_p3
#SBATCH -t 01:00:00
#SBATCH --array=0-5           # Array Index Mapping: 0=mllw, 1=mhhw, 2=mhw, 3=mlw, 4=mtl, 5=dtl, 0-5 means all
#SBATCH --cpus-per-task=4     # 4 CPUs
#SBATCH --mem=64G             # 32GB RAM
#SBATCH -o logs/svu_%A_%a.out
#SBATCH -e logs/svu_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov

set -euo pipefail
# Change this path to your actual environment
source /home/mojganr/miniconda3/bin/activate TAD

# --- IMPORTANT: CPU Threads ---
# This ensures numpy uses all 4 CPUs for matrix operations
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- Configuration ---
RUNID="Pac_SAL"
PRE="/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/pre_process"
OUT="/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/out_nc"

mkdir -p $OUT
mkdir -p logs

echo "=== Running SVU Exact Task ${SLURM_ARRAY_TASK_ID} for RUNID: $RUNID ==="

# Run Python with the Array Index
python3 svu_phase3_pacific.py \
    --runid "${RUNID}" \
    --path_pre "${PRE}" \
    --path_out "${OUT}" \
    --array_index ${SLURM_ARRAY_TASK_ID}
    
echo "Task ${SLURM_ARRAY_TASK_ID} Complete"
