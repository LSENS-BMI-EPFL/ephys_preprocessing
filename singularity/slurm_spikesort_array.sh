#!/bin/bash
#SBATCH --job-name=ephys-array
#SBATCH --array=1-40%10            # 40 tasks, max 8 concurrent (2 H100 nodes)
#SBATCH --output=logs/spikesort_%A_%a.out  # %A=job ID, %a=array task ID
#SBATCH --error=logs/spikesort_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --partition=h100
#SBATCH --gres=gpu:1              # Each task gets 1 GPU
#SBATCH --cpus-per-task=16         # 32 cores / 4 GPUs = 8 cores per job
#SBATCH --mem=32G                 # Conservative: 16GB per job

# ============================================================================
# SLURM job array script for parallel spike sorting (H100 GPUs)
# Each array task processes one session from inputs.txt
# ============================================================================

# Load required modules (adjust for your HPC environment)
# module load singularity
# module load cuda/12.1

# Set environment variables
export EPHYS_LOG_DIR=/scratch/lebert/ephys/log
export CUDA_VISIBLE_DEVICES=0

# Prevent KMeans/MKL deadlocks during Kilosort4 drift correction
# (scikit-learn + CUDA can hang on HPC clusters with multiple threads)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Define paths
SIF_IMAGE="singularity/ephys-pipeline.sif"
CONFIG_DIR="config"
LOG_DIR="logs"
INPUT_FILE="${CONFIG_DIR}/inputs.txt"
CODE_DIR="/home/lebert/code/ephys_preprocessing"   # adjust to HPC clone location

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
echo "GPUs: ${SLURM_GPUS}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Session: $SESSION"
echo "Start Time: $(date)"
echo "=================================="

# Run spike sorting pipeline with GPU support for single session
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
from scripts.preprocess_spikesort_si import main
import yaml

# Load config to get raw_data_path
config_path = Path('/mnt/config/preprocess_config_si_hpc_array.yaml')
with open(config_path) as f:
    config = yaml.safe_load(f)
data_path = Path(config['raw_data_path'])

# Process this session
input_path = '${SESSION}'
print(f'Processing: {data_path / input_path}')
main(data_path / input_path, config_path)
"

# Capture exit status
EXIT_STATUS=$?

# Print completion info
echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: ${EXIT_STATUS}"
echo "=================================="

exit ${EXIT_STATUS}
