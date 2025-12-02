#!/bin/bash
#SBATCH -A vdatum
#SBATCH -J coef_dtl_matlabAligned
#SBATCH -t 1:00:00
#SBATCH --array=0-399
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov

set -euo pipefail

# Activate env
source /home/mojganr/miniconda3/bin/activate TAD

# Keep math libs single-threaded inside each process
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Inputs
HH_NC="/work2/noaa/vdatum/mojganr/work_adcirc/Alaska/a76_9rivers_CF4_TAD_hh.nc"
LL_NC="/work2/noaa/vdatum/mojganr/work_adcirc/Alaska/a76_9rivers_CF4_TAD_ll.nc"
GAGE_CSV="/work2/noaa/vdatum/mojganr/work_adcirc/Alaska/gageinfo_all.csv"

OUTDIR="./out_dtl_matlabAligned_multi"
OUT_NC="coef_full_dtl_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"

mkdir -p "$OUTDIR" logs

# Domain/array config
NCHUNKS="${SLURM_ARRAY_TASK_COUNT:-400}"   # total shards
WORKERS="${SLURM_CPUS_PER_TASK}"
BATCH=20000
MAX_GAP_DAYS=3.0
MAX_STATIONS=301

echo "[$(date '+%F %T')] datum=dtl shard=${SLURM_ARRAY_TASK_ID}/${NCHUNKS} workers=${WORKERS}"

srun -c "${WORKERS}" --cpu-bind=cores \
  python3 S3_Cor_Coe_mp_3_multistations_dtl_mtl.py \
    --hh "$HH_NC" \
    --ll "$LL_NC" \
    --gage "$GAGE_CSV" \
    --datum dtl \
    --outdir "$OUTDIR" \
    --out-nc "$OUT_NC" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --nchunks  "${NCHUNKS}" \
    --batch    "${BATCH}" \
    --workers  "${WORKERS}" \
    --gap-max-days "${MAX_GAP_DAYS}" \
    --max-stations "${MAX_STATIONS}"

echo "[$(date '+%F %T')] DONE datum=dtl shard=${SLURM_ARRAY_TASK_ID}"
