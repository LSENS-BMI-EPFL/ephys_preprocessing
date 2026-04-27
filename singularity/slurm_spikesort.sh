#!/bin/bash
#SBATCH --job-name=ephys-spikesort
#SBATCH --output=logs/spikesort_%j.out
#SBATCH --error=logs/spikesort_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=h100          # Adjust to your cluster's GPU partition name
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# ============================================================================
# SLURM job script for spike sorting pipeline (requires GPU for Kilosort)
# ============================================================================

# Load required modules (adjust for your HPC environment)
# module load singularity
# module load cuda/12.1

# Set environment variables
export EPHYS_LOG_DIR=/opt/ephys/log
export CUDA_VISIBLE_DEVICES=0

# Define paths
SIF_IMAGE="singularity/ephys-pipeline.sif"
CONFIG_DIR="config"
LOG_DIR="logs"
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

# Print job info
echo "=================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Node: ${SLURM_NODELIST}"
echo "Start Time: $(date)"
echo "=================================="

# Run spike sorting pipeline with GPU support
singularity exec \
  --nv \
  --bind ${BIND_DATA} \
  --bind ${BIND_OUTPUT} \
  --bind ${BIND_MICE} \
  --bind ${BIND_CONFIG} \
  --bind ${BIND_LOGS} \
  --bind ${BIND_CODE} \
  ${SIF_IMAGE} \
  python3.11 /opt/ephys/scripts/preprocess_spikesort_si.py \
    --input-list /mnt/config/inputs.txt \
    --config /mnt/config/preprocess_config_si_hpc.yaml

# Print completion info
echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: $?"
echo "=================================="
