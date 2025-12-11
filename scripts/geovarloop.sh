#!/bin/bash
# Simple sequential runner for geovar_all_chroms.py

# activate your env
module load conda
eval "$(mamba shell hook --shell bash)"
mamba activate /quobyte/bmhenngrp/conda_envs/python_geovar

# path to your script
PY_SCRIPT="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/scripts/geovar_all_chroms.py"

# make logs folder
LOG_DIR="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/scripts/logs"
mkdir -p "$LOG_DIR"

#list of chromosomes (edit if needed)
for CHR in {1..22}
do
    echo "[INFO] Starting chr${CHR} at $(date)"
    python "$PY_SCRIPT" --chr "$CHR" &> "${LOG_DIR}/geovar_chr${CHR}.log"
    echo "[INFO] Finished chr${CHR} at $(date)"
done

echo "all is done dude!"
