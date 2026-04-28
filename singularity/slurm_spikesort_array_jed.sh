#!/bin/bash
#SBATCH --job-name=ephys-array-cpu
#SBATCH --array=1-40%10            # 40 tasks, max 10 concurrent
#SBATCH --output=logs/spikesort_%A_%a.out
#SBATCH --error=logs/spikesort_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --qos=parallel

# ============================================================================

# SLURM job array script for parallel spike sorting (CPU version)

# Each array task processes one session from inputs.txt

# ============================================================================

# Load required modules (adjust if needed)

# module load singularity

# Set environment variables

export EPHYS_LOG_DIR=/home/bisi/logs

# Prevent thread oversubscription / deadlocks

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Define paths

SIF_IMAGE="singularity/ephys-pipeline.sif"
CONFIG_DIR="config"
LOG_DIR="logs"
INPUT_FILE="${CONFIG_DIR}/inputs.txt"
CODE_DIR="/home/bisi/code/ephys_preprocessing"

# Singularity bind paths

BIND_DATA="/scratch/bisi/data:/scratch/bisi/data:ro"
BIND_OUTPUT="/scratch/bisi/ephys_output:/scratch/bisi/ephys_output"
BIND_MICE="/home/bisi/mice_info:/home/bisi/mice_info:ro"
BIND_CONFIG="${CONFIG_DIR}:/mnt/config:ro"
BIND_LOGS="${LOG_DIR}:/opt/ephys/log"
BIND_CODE="${CODE_DIR}/ephys_preprocessing:/opt/ephys/ephys_preprocessing,${CODE_DIR}/scripts:/opt/ephys/scripts"

# Create log directory if it doesn't exist

mkdir -p ${LOG_DIR}

# Get session from inputs.txt

SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE")

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

singularity exec
--bind ${BIND_DATA}
--bind ${BIND_OUTPUT}
--bind ${BIND_MICE}
--bind ${BIND_CONFIG}
--bind ${BIND_LOGS}
--bind ${BIND_CODE}
${SIF_IMAGE}
python3.11 -c "
import sys
sys.path.insert(0, '/opt/ephys')
from pathlib import Path
from scripts.preprocess_spikesort_si import main
import yaml

config_path = Path('/mnt/config/preprocess_config_si_hpc_array_axel.yaml')
with open(config_path) as f:
config = yaml.safe_load(f)
data_path = Path(config['raw_data_path'])

input_path = '${SESSION}'
print(f'Processing: {data_path / input_path}')
main(data_path / input_path, config_path)
"

EXIT_STATUS=$?

echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: ${EXIT_STATUS}"
echo "=================================="

exit ${EXIT_STATUS}
