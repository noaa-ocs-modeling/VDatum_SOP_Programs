#!/bin/bash
#SBATCH -J svu_p2
#SBATCH -t 4:00:00
#SBATCH --array=0-398
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o logs/svu_p2_%A_%a.out
#SBATCH -e logs/svu_p2_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov
set -euo pipefail
source /home/mojganr/miniconda3/bin/activate TAD

RUNID="Pac_SAL"
DATUMS=("mllw" "mhhw" "mhw" "mlw" "mtl" "dtl")
#DATUMS=("mllw")
PRE="/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/pre_process/"
CHUNKS=399

cd /work2/noaa/vdatum/mojganr/work_adcirc/SVU/S4_svu/Pacific/

echo "=== Running Phase 2 chunks for all datums ==="

for VAR in "${DATUMS[@]}"; do
    python3 svu_phase2_pacific.py \
        --runid $RUNID \
        --datum $VAR \
        --chunk-id $SLURM_ARRAY_TASK_ID \
        --nchunks $CHUNKS \
        --path_pre $PRE
done
    
    

