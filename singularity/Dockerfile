# Dockerfile for ephys preprocessing pipeline
# Supports both Docker and Singularity/Apptainer (via conversion)

# ============================================================================
# Stage 1: Builder stage - Install dependencies
# ============================================================================
FROM nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04 AS builder

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

# FIX: Pin PyTorch to CUDA 12.1 BEFORE installing other deps via uv.
# Without this explicit pin, uv/pip resolves the latest PyTorch wheel
# which is built against CUDA 13.x — incompatible with kuma's driver
# which supports up to CUDA 12.6 (driver version 12060).
RUN python3.11 -m pip install --no-cache-dir \
    torch \
    --index-url https://download.pytorch.org/whl/cu121

# Install remaining Python dependencies using uv.
# PyTorch is already installed above so it will not be overridden.
RUN python3.11 -m uv pip install --system --no-cache .

# Verify PyTorch CUDA version at build time — fails the build if wrong.
RUN python3.11 -c "\
import torch; \
v = torch.version.cuda; \
assert v.startswith('12.1'), f'Expected CUDA 12.1, got {v}. Check pyproject.toml for a conflicting torch pin.'; \
print(f'PyTorch CUDA version OK: {v}')"

# ============================================================================
# Stage 2: Runtime stage - Minimal production image
# ============================================================================
# Keep cudnn8-runtime for a smaller image. Base must stay 12.1.x to match
# the PyTorch wheel installed above and remain under kuma's driver cap (12.6).
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
ENV EPHYS_LOG_DIR=/opt/ephys/log

# CUDA environment variables (for PyTorch/Kilosort)
# NOTE: CUDA_VISIBLE_DEVICES is intentionally NOT set here.
# SLURM controls GPU assignment via --gres=gpu:1. Setting it in the image
# can conflict with SLURM's assignment on multi-GPU nodes.
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
LABEL version="0.2.0"

# ============================================================================
# Build, push, and verify instructions:
# ============================================================================
# 1. Build:
#      docker build -t your_dockerhub_username/ephys-pipeline:cu121-v0.2.0 .
#      docker tag your_dockerhub_username/ephys-pipeline:cu121-v0.2.0 \
#                 your_dockerhub_username/ephys-pipeline:latest
#
# 2. Verify GPU locally before pushing:
#      docker run --rm --runtime=nvidia your_dockerhub_username/ephys-pipeline:latest \
#        python3.11 -c "import torch; print(torch.cuda.is_available(), torch.version.cuda)"
#      -> Expected: True 12.1
#
# 3. Push both tags:
#      docker push your_dockerhub_username/ephys-pipeline:cu121-v0.2.0
#      docker push your_dockerhub_username/ephys-pipeline:latest
#
# 4. Pull on kuma (use YOUR image, not alpine):
#      apptainer pull ephys-pipeline_latest.sif \
#        docker://your_dockerhub_username/ephys-pipeline:latest
#
# 5. Verify on kuma before submitting jobs:
#      apptainer exec --nv ephys-pipeline_latest.sif \
#        python3.11 -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
#      -> Expected: True  NVIDIA H100 ...
