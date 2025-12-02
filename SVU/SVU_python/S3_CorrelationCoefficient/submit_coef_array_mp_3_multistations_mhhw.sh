#!/bin/bash
#SBATCH -A vdatum
#SBATCH -J coef_mhhw_matlabAligned
#SBATCH -t 2:00:00
#SBATCH --array=0-399
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov

set -euo pipefail

# Activate environment
source /home/mojganr/miniconda3/bin/activate TAD

# Keep math libs single-threaded inside each Python worker
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Inputs
#HH_NC="/work2/noaa/vdatum/mojganr/work_adcirc/Alaska/a76_9rivers_CF4_TAD_hh.nc"
#GAGE_CSV="/work2/noaa/vdatum/mojganr/work_adcirc/Alaska/gageinfo_all.csv"
HH_NC="/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_08182025_TAD/TADhh.nc"
GAGE_CSV="/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/gageinfo.csv"



# Output (per shard)
OUTDIR="./out_mhhw_matlabAligned_multi"
OUT_NC="coef_full_mhhw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"
mkdir -p "$OUTDIR" logs

# Sharding / workers
NCHUNKS="${SLURM_ARRAY_TASK_COUNT:-400}"
WORKERS="${SLURM_CPUS_PER_TASK}"
BATCH=20000
MAX_GAP_DAYS=3.0
#MAX_STATIONS=301  #Alaska
MAX_STATIONS=60    
echo "[$(date '+%F %T')] shard=${SLURM_ARRAY_TASK_ID}/${NCHUNKS} workers=${WORKERS}"

srun -c "${WORKERS}" --cpu-bind=cores \
  python3 S3_Cor_Coe_mp_3_multistations_mhhw.py \
    --hh  "$HH_NC" \
    --gage "$GAGE_CSV" \
    --outdir "$OUTDIR" \
    --out-nc "$OUT_NC" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --nchunks  "${NCHUNKS}" \
    --batch    "${BATCH}" \
    --workers  "${WORKERS}" \
    --gap-max-days "${MAX_GAP_DAYS}"\
    --max-stations "${MAX_STATIONS}"

echo "[$(date '+%F %T')] DONE shard=${SLURM_ARRAY_TASK_ID}"
