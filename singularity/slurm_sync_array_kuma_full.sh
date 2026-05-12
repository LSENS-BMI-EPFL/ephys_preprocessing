#!/bin/bash
#SBATCH --job-name=ephys-array
#SBATCH --array=1-3%3           # 40 tasks, max 8 concurrent (2 H100 nodes)
#SBATCH --output=logs/sync_%A_%a.out
#SBATCH --error=logs/sync_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=h100
#SBATCH --gres=gpu:1              # Each task gets 1 GPU
#SBATCH --cpus-per-task=16         # 64 cores / 4 GPUs = 16 cores per job
#SBATCH --mem=48G

#SBATCH --mail-user=axel.bisi@epfl.ch
#SBATCH --mail-type=END
# ============================================================================
# SLURM job array script for parallel sync processing (CPU only)
# Each array task processes one session from inputs.txt
# Runs TPrime, Cwaves, and waveform metrics
# ============================================================================

# Load required modules (adjust for your HPC environment)
# module load singularity

# Set environment variables
# EPHYS_LOG_DIR=/scratch/bisi/ephys/log
export EPHYS_LOG_DIR=/home/bisi/logs

# Define paths
SIF_IMAGE="ephys-pipeline_latest.sif" # currently in folder singularity
CONFIG_DIR="/home/bisi/code/ephys_preprocessing/config"
LOG_DIR="/home/bisi/logs"
INPUT_FILE="${CONFIG_DIR}/inputs_axel_test.txt" # note
CODE_DIR="/home/bisi/code/ephys_preprocessing"

# Singularity bind paths
BIND_DATA="/scratch/bisi/data:/scratch/bisi/data"
BIND_OUTPUT="/scratch/bisi/data:/scratch/bisi/data"
BIND_MICE="/home/bisi/mice_info:/home/bisi/mice_info:ro"
BIND_CONFIG="${CONFIG_DIR}:/mnt/config:ro"
BIND_LOGS="${LOG_DIR}:/opt/ephys/log"
BIND_CODE="${CODE_DIR}/ephys_preprocessing:/opt/ephys/ephys_preprocessing,${CODE_DIR}/scripts:/opt/ephys/scripts"

# Create log directory if it doesn't exist
mkdir -p ${LOG_DIR}

# Get the Nth line from inputs.txt (where N = SLURM_ARRAY_TASK_ID)
#SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE")
SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE" | tr -d '\n\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | sed 's:/*$::')

# Skip if empty or comment
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

# Run sync processing pipeline (NO --nv flag, no GPU)
# Transforms input path to find catgt folder automatically
singularity exec \
  --nv \
  --bind ${BIND_DATA} \
  --bind ${BIND_OUTPUT} \
  --bind ${BIND_MICE} \
  --bind ${BIND_CONFIG} \
  --bind ${BIND_LOGS} \
  --bind ${BIND_CODE} \
  ${SIF_IMAGE} \
  bash -c "
python -m pip install --user ibllib ibl-neuropixel
"

singularity exec \
  --nv \
  --bind ${BIND_DATA} \
  --bind ${BIND_OUTPUT} \
  --bind ${BIND_MICE} \
  --bind ${BIND_CONFIG} \
  --bind ${BIND_LOGS} \
  --bind ${BIND_CODE} \
  ${SIF_IMAGE} \
  python3.11 -c "
import sys
sys.path.insert(0, '/opt/ephys')
from pathlib import Path
from scripts.preprocess_sync_si import main, transform_input_to_catgt_path
import yaml


# Load config
config_path = Path('/mnt/config/preprocess_config_si_hpc_array.yaml')
with open(config_path) as f:
    config = yaml.safe_load(f)

# Use output_path as the data root (where catgt folders are)
data_root = Path(config['output_path'])

# Transform input path to find catgt folder
input_path = '${SESSION}'
try:
    catgt_path = transform_input_to_catgt_path(input_path, data_root)
    print(f'Processing: {catgt_path}')
    print(f'  (transformed from: {input_path})')
    main(catgt_path, config)
except FileNotFoundError as e:
    print(f'ERROR: {e}')
    print(f'Skipping session: {input_path}')
    sys.exit(1)
"

# Capture exit status
EXIT_STATUS=$?

# Print completion info
echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: ${EXIT_STATUS}"
echo "=================================="

exit ${EXIT_STATUS}
