#!/bin/bash
#SBATCH -A vdatum
#SBATCH -J coef_all_datums_matlabAligned_multi
#SBATCH -t 2:00:00
#SBATCH --array=0-399
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mojgan.rostaminia@noaa.gov

set -euo pipefail
# --------- go to code directory (IMPORTANT) ----------
CDIR="/work2/noaa/vdatum/mojganr/work_adcirc/SVU/S3_CC"
cd "$CDIR"

source /home/mojganr/miniconda3/bin/activate  TAD

# Keep math libs single-threaded inside each Python worker
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# ---------------- INPUT FILES ----------------
HH_NC="/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_08182025_TAD/TADhh.nc"
LL_NC="/work2/noaa/vdatum/mojganr/work_adcirc/TAD/Pacific_TAD/PA_6sec_SAL_08182025_TAD/TADll.nc"
GAGE_CSV="/work2/noaa/vdatum/mojganr/work_adcirc/Pacific/gageinfo.csv"

# ---------------- PARALLEL CONFIG ----------------
NCHUNKS="${SLURM_ARRAY_TASK_COUNT:-400}"   # total shards across array
WORKERS="${SLURM_CPUS_PER_TASK}"
BATCH=20000
MAX_GAP_DAYS=3.0
MAX_STATIONS=60

mkdir -p logs

# ---------------- PYTHON SCRIPTS (EDIT NAMES IF NEEDED) ----------------
# These are the *Python files* you already have for each datum.
# Rename these variables if your filenames are slightly different.
SCRIPT_MLW="S3_Cor_Coe_mp_3_multistations_mlw.py"   
SCRIPT_MLLW="S3_Cor_Coe_mp_3_multistations_mllw.py"
SCRIPT_MHW="S3_Cor_Coe_mp_3_multistations_mhw.py"
SCRIPT_MHHW="S3_Cor_Coe_mp_3_multistations_mhhw.py"
SCRIPT_MTL_DTL="S3_Cor_Coe_mp_3_multistations_dtl_mtl.py"   # handles both mtl & dtl

# ---------------- DATUM LOOP ----------------
DATUMS=(mlw mllw mhw mhhw mtl dtl)

for D in "${DATUMS[@]}"; do
  echo "[$(date '+%F %T')] START datum=${D} shard=${SLURM_ARRAY_TASK_ID}/${NCHUNKS}"

  case "$D" in
    mlw)
      OUTDIR="./out_mlw_matlabAligned_multi"
      OUT_NC="coef_full_mlw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"
      mkdir -p "$OUTDIR"

      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MLW" \
          --ll "$LL_NC" \
          --gage "$GAGE_CSV" \
          --outdir "$OUTDIR" \
          --out-nc "$OUT_NC" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;

    mllw)
      OUTDIR="./out_mllw_matlabAligned_multi"
      OUT_NC="coef_full_mllw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"
      mkdir -p "$OUTDIR"
      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MLLW" \
          --ll "$LL_NC" \
          --gage "$GAGE_CSV" \
          --outdir "$OUTDIR" \
          --out-nc "$OUT_NC" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;

    mhw)
      OUTDIR="./out_mhw_matlabAligned_multi"
      OUT_NC="coef_full_mhw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"
      mkdir -p "$OUTDIR"
      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MHW" \
          --hh "$HH_NC" \
          --gage "$GAGE_CSV" \
          --outdir "$OUTDIR" \
          --out-nc "$OUT_NC" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;

    mhhw)
      OUTDIR="./out_mhhw_matlabAligned_multi"
      OUT_NC="coef_full_mhhw_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc"
      mkdir -p "$OUTDIR"
      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MHHW" \
          --hh "$HH_NC" \
          --gage "$GAGE_CSV" \
          --outdir "$OUTDIR" \
          --out-nc "$OUT_NC" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;

    mtl)
      OUTDIR="./out_mtl_matlabAligned_multi"
      mkdir -p "$OUTDIR"
      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MTL_DTL" \
          --hh "$HH_NC" \
          --ll "$LL_NC" \
          --gage "$GAGE_CSV" \
          --datum mtl \
          --outdir "$OUTDIR" \
          --out-nc "coef_full_mtl_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;

    dtl)
      OUTDIR="./out_dtl_matlabAligned_multi"
      mkdir -p "$OUTDIR"
      srun -c "${WORKERS}" --cpu-bind=cores \
        python3 "$SCRIPT_MTL_DTL" \
          --hh "$HH_NC" \
          --ll "$LL_NC" \
          --gage "$GAGE_CSV" \
          --datum dtl \
          --outdir "$OUTDIR" \
          --out-nc "coef_full_dtl_matlabAligned_multi_SLICE_${SLURM_ARRAY_TASK_ID}.nc" \
          --task-id "${SLURM_ARRAY_TASK_ID}" \
          --nchunks "${NCHUNKS}" \
          --batch "${BATCH}" \
          --workers "${WORKERS}" \
          --gap-max-days "${MAX_GAP_DAYS}" \
          --max-stations "${MAX_STATIONS}"
      ;;
  esac

  echo "[$(date '+%F %T')] DONE  datum=${D} shard=${SLURM_ARRAY_TASK_ID}"
done

