#!/bin/bash
#SBATCH --job-name=ephys-array
#SBATCH --array=1-8%3
#SBATCH --output=logs/ephys_%A_%a.out
#SBATCH --error=logs/ephys_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=h100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --mail-user=axel.bisi@epfl.ch
#SBATCH --mail-type=END
# ============================================================================
# SLURM job array — unified ephys pipeline (spikesort + sync) per session
# Each task reads one line from inputs_axel.txt and runs preprocess_si.py
# ============================================================================

export EPHYS_LOG_DIR=/home/bisi/logs
export CUDA_VISIBLE_DEVICES=0

# Prevent KMeans/MKL deadlocks during Kilosort4 drift correction
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

SIF_IMAGE="ephys-pipeline_latest.sif"
CONFIG_DIR="/home/bisi/code/ephys_preprocessing/config"
CONFIG_FILE="${CONFIG_DIR}/preprocess_config_si_hpc_array_axel_full.yaml"
LOG_DIR="/home/bisi/logs"
INPUT_FILE="${CONFIG_DIR}/inputs_axel_test.txt"
CODE_DIR="/home/bisi/code/ephys_preprocessing"

BIND_DATA="/scratch/bisi/data:/scratch/bisi/data"
BIND_MICE="/home/bisi/mice_info:/home/bisi/mice_info:ro"
BIND_CONFIG="${CONFIG_DIR}:/mnt/config:ro"
BIND_LOGS="${LOG_DIR}:/opt/ephys/log"
BIND_CODE="${CODE_DIR}/ephys_preprocessing:/opt/ephys/ephys_preprocessing,${CODE_DIR}/scripts:/opt/ephys/scripts"

mkdir -p "${LOG_DIR}"

# Read and sanitise the session path for this array task
SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INPUT_FILE}" \
          | tr -d '\n\r' \
          | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' \
          | sed 's:/*$::')

if [ -z "$SESSION" ] || [[ "$SESSION" =~ ^[[:space:]]*# ]]; then
    echo "No valid session for array task ${SLURM_ARRAY_TASK_ID}"
    echo "Line content: '${SESSION}'"
    exit 0
fi

echo "=================================="
echo "Array Job ID : ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Node         : ${SLURM_NODELIST}"
echo "CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "Session      : ${SESSION}"
echo "Start Time   : $(date)"
echo "=================================="

singularity exec \
    --nv \
    --bind ${BIND_DATA} \
    --bind ${BIND_MICE} \
    --bind ${BIND_CONFIG} \
    --bind ${BIND_LOGS} \
    --bind ${BIND_CODE} \
    ${SIF_IMAGE} \
    python3.11 -c "
import sys
sys.path.insert(0, '/opt/ephys')
from pathlib import Path
from scripts.preprocess_si import main
import yaml

config_path = Path('/mnt/config/preprocess_config_si_hpc_array_axel_full.yaml')
session = Path('''${SESSION}''')

if not session.is_absolute():
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    session = Path(cfg['raw_data_path']) / session

print(f'Processing: {session}')
main(session, config_path)
"

EXIT_STATUS=$?

echo "=================================="
echo "End Time   : $(date)"
echo "Exit Status: ${EXIT_STATUS}"
echo "=================================="

exit ${EXIT_STATUS}
