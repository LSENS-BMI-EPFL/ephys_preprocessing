# Dockerfile for ephys preprocessing pipeline
# Supports both Docker and Singularity (via conversion)

# ============================================================================
# Stage 1: Builder stage - Install dependencies
# ============================================================================
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04 AS builder

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    python3-pip \
    git \
    wget \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.11 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

# Upgrade pip and install uv for faster dependency resolution
RUN python3.11 -m pip install --no-cache-dir --upgrade pip && \
    python3.11 -m pip install --no-cache-dir uv

# Create working directory
WORKDIR /opt/ephys

# Copy only dependency files first (for layer caching)
COPY pyproject.toml ./

# Install Python dependencies using uv
RUN python3.11 -m uv pip install --system --no-cache .

# ============================================================================
# Stage 2: Runtime stage - Minimal production image
# ============================================================================
FROM nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Enable 32-bit architecture support for CatGT and other tools
RUN dpkg --add-architecture i386

# Install only runtime dependencies (including 32-bit libraries for CatGT)
RUN apt-get update && apt-get install -y \
    python3.11 \
    libqt5core5a \
    libqt5gui5 \
    libqt5widgets5 \
    libxcb-xinerama0 \
    libgomp1 \
    libc6:i386 \
    libstdc++6:i386 \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.11 as default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

# Copy Python environment from builder
COPY --from=builder /usr/local/lib/python3.11 /usr/local/lib/python3.11
COPY --from=builder /usr/local/bin /usr/local/bin

# Set up working directory
WORKDIR /opt/ephys

# ============================================================================
# Download and install external Linux binaries from SpikeGLX
# ============================================================================
RUN mkdir -p /opt/tools && cd /opt/tools && \
    # Download CatGT
    wget -q https://billkarsh.github.io/SpikeGLX/Support/CatGTLnxApp.zip && \
    unzip -q CatGTLnxApp.zip && \
    rm CatGTLnxApp.zip && \
    cd CatGT-linux && \
    chmod +x install.sh && ./install.sh && \
    # Download TPrime
    cd /opt/tools && \
    wget -q https://billkarsh.github.io/SpikeGLX/Support/TPrimeLnxApp.zip && \
    unzip -q TPrimeLnxApp.zip && \
    rm TPrimeLnxApp.zip && \
    cd TPrime-linux && \
    chmod +x install.sh && ./install.sh && \
    # Download C_Waves
    cd /opt/tools && \
    wget -q https://billkarsh.github.io/SpikeGLX/Support/C_WavesLnxApp.zip && \
    unzip -q C_WavesLnxApp.zip && \
    rm C_WavesLnxApp.zip && \
    cd C_Waves-linux && \
    chmod +x install.sh && ./install.sh && \
    # Download OverStrike
    cd /opt/tools && \
    wget -q https://billkarsh.github.io/SpikeGLX/Support/OverStrikeLnxApp.zip && \
    unzip -q OverStrikeLnxApp.zip && \
    rm OverStrikeLnxApp.zip && \
    cd OverStrike-linux && \
    chmod +x install.sh && ./install.sh

# ============================================================================
# Copy Python package
# ============================================================================
COPY ephys_preprocessing /opt/ephys/ephys_preprocessing
COPY scripts /opt/ephys/scripts

# Create log directory
RUN mkdir -p /opt/ephys/log

# ============================================================================
# Environment variables
# ============================================================================
ENV PYTHONPATH=/opt/ephys
ENV PATH=/opt/tools/CatGT-linux:/opt/tools/TPrime-linux:/opt/tools/OverStrike-linux:/opt/tools/C_Waves-linux:$PATH

# Set log directory environment variable
ENV EPHYS_LOG_DIR=/opt/ephys/log

# CUDA environment variables (for PyTorch/Kilosort)
ENV CUDA_VISIBLE_DEVICES=0
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

# ============================================================================
# Entrypoint and command
# ============================================================================
ENTRYPOINT ["python3.11", "/opt/ephys/scripts/preprocess_spikesort_si.py"]
CMD ["--help"]

# ============================================================================
# Labels for metadata
# ============================================================================
LABEL maintainer="Jules Lebert"
LABEL description="Electrophysiology preprocessing pipeline with SpikeInterface, Kilosort4, and Bombcell"
LABEL version="0.1.0"

# ============================================================================
# Build instructions:
# ============================================================================
# Build the image (tools are downloaded automatically):
# docker build -t ephys-pipeline:latest .
#
# Or use docker compose:
# docker compose build
#
# Test the build:
# docker run --rm --runtime=nvidia ephys-pipeline:latest --help
#
# Run with data:
# docker run --runtime=nvidia \
#   -v /mnt/lsens/data:/mnt/lsens/data:ro \
#   -v /home/lebert/lsens_srv/analysis/Jules_Lebert/data:/home/lebert/lsens_srv/analysis/Jules_Lebert/data \
#   -v $(pwd)/config:/mnt/config:ro \
#   ephys-pipeline:latest \
#   --input-list /mnt/config/inputs.txt \
#   --config /mnt/config/preprocess_config_si_docker.yaml
#
# Or use docker compose:
# docker compose run --rm ephys-spikesort
