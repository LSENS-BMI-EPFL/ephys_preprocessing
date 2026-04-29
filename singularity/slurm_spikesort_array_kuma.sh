#!/bin/bash
#SBATCH --job-name=ephys-array
#SBATCH --array=1-3%10            # 40 tasks, max 8 concurrent (2 H100 nodes)
#SBATCH --output=logs/spikesort_%A_%a.out
#SBATCH --error=logs/spikesort_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=h100
#SBATCH --gres=gpu:1              # Each task gets 1 GPU
#SBATCH --cpus-per-task=16         # 32 cores / 4 GPUs = 8 cores per job
#SBATCH --mem=32G
#SBATCH --qos=long
# ============================================================================

# SLURM job array script for parallel spike sorting (CPU version)

# Each array task processes one session from inputs.txt

# ============================================================================

# Load required modules (adjust if needed)
# module load singularity
# module load cuda/12.1
# module load singularity

# Set environment variables
# Set environment variables
export EPHYS_LOG_DIR=/home/bisi/logs
export CUDA_VISIBLE_DEVICES=0

# Prevent KMeans/MKL deadlocks during Kilosort4 drift correction
# (scikit-learn + CUDA can hang on HPC clusters with multiple threads)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Define paths

SIF_IMAGE="ephys-pipeline_latest.sif" # currently in folder singularity
CONFIG_DIR="/home/bisi/code/ephys_preprocessing/config"
LOG_DIR="/home/bisi/logs"
INPUT_FILE="${CONFIG_DIR}/inputs_axel_test.txt" # note
CODE_DIR="/home/bisi/code/ephys_preprocessing"

# Singularity bind paths

BIND_DATA="/scratch/bisi/data:/scratch/bisi/data" #removed-read only
BIND_OUTPUT="/scratch/bisi/data:/scratch/bisi/data"
BIND_MICE="/home/bisi/mice_info:/home/bisi/mice_info:ro"
BIND_CONFIG="${CONFIG_DIR}:/mnt/config:ro"
BIND_LOGS="${LOG_DIR}:/opt/ephys/log"
BIND_CODE="${CODE_DIR}/ephys_preprocessing:/opt/ephys/ephys_preprocessing,${CODE_DIR}/scripts:/opt/ephys/scripts"

# Create log directory if it doesn't exist

mkdir -p ${LOG_DIR}

# Get session from inputs.txt
echo "Reading session for array task ${SLURM_ARRAY_TASK_ID} from $(realpath ${INPUT_FILE})"
#SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE")
#SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE" | xargs)
SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE" | tr -d '\n\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sed 's:/*$::')

# Skip invalid lines

if [ -z "$SESSION" ] || [[ "$SESSION" =~ ^[[:space:]]*# ]]; then
echo "No valid session for array task ${SLURM_ARRAY_TASK_ID}"
echo "Line content: '$SESSION'"
exit 0
fi

# Print job info

echo "=================================="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Session: $SESSION"
echo "Start Time: $(date)"
echo "=================================="

# Run pipeline (CPU only, no --nv)

singularity exec \
--bind ${BIND_DATA} \
--bind ${BIND_OUTPUT} \
--bind ${BIND_MICE} \
--bind ${BIND_CONFIG} \
--bind ${BIND_LOGS} \
--bind ${BIND_CODE} \
${SIF_IMAGE} \
python3.11 -c "
import sys
import os
sys.path.insert(0, '/opt/ephys')
from pathlib import Path
from scripts.preprocess_spikesort_si import main
import yaml

config_path = Path('/mnt/config/preprocess_config_si_hpc_array_axel.yaml')
with open(config_path) as f:
    config = yaml.safe_load(f)

input_path = Path('''${SESSION}''')

# Use full path if it's absolute, otherwise join with data_path
if input_path.is_absolute():
    full_path = input_path
else:
    data_path = Path(config['raw_data_path'])
    full_path = data_path / input_path

print(f'Processing: {full_path}')
main(full_path, config_path)
"

EXIT_STATUS=$?

echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: ${EXIT_STATUS}"
echo "=================================="

exit ${EXIT_STATUS}
