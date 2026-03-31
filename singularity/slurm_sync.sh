#!/bin/bash
#SBATCH --job-name=ephys-sync
#SBATCH --output=logs/sync_%j.out
#SBATCH --error=logs/sync_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=cpu          # Adjust to your cluster's CPU partition name
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# ============================================================================
# SLURM job script for sync processing pipeline (CPU only, no GPU needed)
# Runs TPrime, Cwaves, waveform metrics, and LFP analysis
# ============================================================================

# Load required modules (adjust for your HPC environment)
# module load singularity

# Set environment variables
export EPHYS_LOG_DIR=/opt/ephys/log

# Define paths
SIF_IMAGE="singularity/ephys-pipeline.sif"
CONFIG_DIR="config"
LOG_DIR="logs"
CODE_DIR="/path/to/ephys_preprocessing"   # adjust to HPC clone location

# Singularity bind paths
BIND_DATA="/home/lebert/lsens_srv/data:/home/lebert/lsens_srv/data:ro"
BIND_OUTPUT="/home/lebert/lsens_srv/analysis/Jules_Lebert/data:/home/lebert/lsens_srv/analysis/Jules_Lebert/data"
BIND_MICE="/home/lebert/lsens_srv/analysis/Jules_Lebert/mice_info:/home/lebert/lsens_srv/analysis/Jules_Lebert/mice_info:ro"
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

# Run sync processing pipeline (NO --nv flag, no GPU)
# Uses same inputs.txt as spikesort - paths are transformed automatically to find catgt folders
singularity exec \
  --bind ${BIND_DATA} \
  --bind ${BIND_OUTPUT} \
  --bind ${BIND_MICE} \
  --bind ${BIND_CONFIG} \
  --bind ${BIND_LOGS} \
  --bind ${BIND_CODE} \
  ${SIF_IMAGE} \
  python3.11 /opt/ephys/scripts/preprocess_sync_si.py \
    --input-list /mnt/config/inputs.txt \
    --config /mnt/config/preprocess_config_si_docker.yaml

# Print completion info
echo "=================================="
echo "End Time: $(date)"
echo "Exit Status: $?"
echo "=================================="
