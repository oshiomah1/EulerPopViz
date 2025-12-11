#!/bin/bash
#SBATCH --job-name=geovar_array
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --array=1-22%11          # 1..22 throttle to 11 at a time
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/scripts/logs/%x-%A_%a.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/scripts/logs/%x-%A_%a.err

set -euo pipefail

# ---- Paths ----
PY_SCRIPT="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/scripts/geovar_all_chroms.py"
ENV_PATH="/quobyte/bmhenngrp/conda_envs/python_geovar"

mkdir -p "$(dirname "$PY_SCRIPT")/logs"

# Map array index -> chromosome label
#CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
CHRS=(2 4 5 8 10 11)

CHR="${CHRS[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "[INFO] $(date) Host: $(hostname)"
echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} -> chr${CHR}"
echo "[INFO] Env: ${ENV_PATH}"
echo "[INFO] Script: ${PY_SCRIPT}"

# ---- Activate env ----
module load conda            
eval "$(mamba shell hook --shell bash)"
mamba activate "${ENV_PATH}"

# Optional: show quick version info once per task
python - <<'PYV'
import sys; print("Python:", sys.version.split()[0])
try:
    import geovar, numpy, pandas, tqdm
    print("geovar OK, numpy/pandas/tqdm OK")
except Exception as e:
    print("Import check failed:", e)
PYV

# ---- Run single chromosome ----
srun -c "${SLURM_CPUS_PER_TASK}" python "${PY_SCRIPT}" --chr "${CHR}"

echo "[INFO] $(date) Done chr${CHR}"
