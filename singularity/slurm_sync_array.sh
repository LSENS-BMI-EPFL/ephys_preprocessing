#!/bin/bash
#SBATCH --job-name=ephys-sync-array
#SBATCH --array=1-30%16           # 30 tasks, max 16 concurrent (CPU only, can run more)
#SBATCH --output=logs/sync_%A_%a.out  # %A=job ID, %a=array task ID
#SBATCH --error=logs/sync_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --partition=cpu           # Adjust to your cluster's CPU partition name
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# ============================================================================
# SLURM job array script for parallel sync processing (CPU only)
# Each array task processes one session from inputs.txt
# Runs TPrime, Cwaves, and waveform metrics
# ============================================================================

# Load required modules (adjust for your HPC environment)
# module load singularity

# Set environment variables
export EPHYS_LOG_DIR=/scratch/lebert/ephys/log

# Define paths
SIF_IMAGE="singularity/ephys-pipeline.sif"
CONFIG_DIR="config"
LOG_DIR="logs"
INPUT_FILE="${CONFIG_DIR}/inputs.txt"
CODE_DIR="/path/to/ephys_preprocessing"   # adjust to HPC clone location

# Singularity bind paths
BIND_DATA="/scratch/lebert/ephys_data:/scratch/lebert/ephys_data:ro"
BIND_OUTPUT="/scratch/lebert/ephys_output:/scratch/lebert/ephys_output"
BIND_MICE="/scratch/lebert/mice_info:/scratch/lebert/mice_info:ro"
BIND_CONFIG="${CONFIG_DIR}:/mnt/config:ro"
BIND_LOGS="${LOG_DIR}:/opt/ephys/log"
BIND_CODE="${CODE_DIR}/ephys_preprocessing:/opt/ephys/ephys_preprocessing,${CODE_DIR}/scripts:/opt/ephys/scripts"

# Create log directory if it doesn't exist
mkdir -p ${LOG_DIR}

# Get the Nth line from inputs.txt (where N = SLURM_ARRAY_TASK_ID)
SESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE")

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
