# Singularity Scripts for HPC

This directory contains scripts for running the ephys preprocessing pipeline on HPC clusters using Singularity.

## Files

- **`ephys-pipeline.sif`** - Singularity image (built from Docker image)
- **`slurm_spikesort.sh`** - SLURM job script for spike sorting (GPU required)
- **`slurm_sync.sh`** - SLURM job script for sync processing (CPU only)
- **`slurm_spikesort_array.sh`** - SLURM job array script for parallel spike sorting
- **`slurm_sync_array.sh`** - SLURM job array script for parallel sync processing

## Quick Start

#### 0. Prerequesites
- Install docker and required tools
- Singularity is on the cluster so no need to install
- 
### 1. Build the Singularity Image

```bash
# From the repo root
docker compose build
singularity build singularity/ephys-pipeline.sif docker-daemon://ephys-pipeline:latest
```
Or, after pushing the Docker image to Docker Hub:
```bash
docker tag source_image:TAG target_image:TAG
docker push username/target_image:TAG
```
then upload an image directly from Dockerhub to Apptainer/Singularity via:
https://scitas-doc.epfl.ch/advanced-guide/singularity-docker/

### 2. Configure for Your HPC Cluster

Edit the SLURM scripts to match your cluster's configuration:

**In all 4 SLURM scripts**, set `CODE_DIR` to the path where the repo is cloned on the cluster:
```bash
CODE_DIR="/path/to/ephys_preprocessing"
```

This bind-mounts the live code into the container at runtime, so you can update the code without rebuilding the `.sif`.

**In `slurm_spikesort.sh` and `slurm_spikesort_array.sh`:**
- Change `#SBATCH --partition=gpu` / `--partition=h100` to your GPU partition name
- Adjust time limits, memory, CPU counts as needed
- Uncomment and adjust `module load` commands if needed

### 3. Prepare Input Files

Create a `config/` directory with:
- **`inputs.txt`** - List of input paths (one per line)
- **`preprocess_config_si_docker.yaml`** - Configuration file

Example `inputs.txt`:
```
JL007/Recording/JL007_20250603_150143/Ephys
JL007/Recording/JL007_20250605_145217/Ephys
```

### 4. Submit Jobs

**Single session:**
```bash
# Submit spike sorting only
sbatch slurm_spikesort.sh

# Submit with dependency (sync runs after spike sorting completes)
JOB1=$(sbatch --parsable slurm_spikesort.sh)
sbatch --dependency=afterok:$JOB1 slurm_sync.sh
```

**Parallel sessions (job arrays):**
```bash
# Submit spike sorting array
sbatch slurm_spikesort_array.sh

# Submit sync array with dependency
JOB1=$(sbatch --parsable slurm_spikesort_array.sh)
sbatch --dependency=afterok:$JOB1 slurm_sync_array.sh
```

## Updating the Code

To deploy code changes to the cluster, just sync the repo and resubmit — no container rebuild needed:

```bash
# On the cluster
git pull

# Or from local machine
rsync -av ephys_preprocessing/ user@cluster:/path/to/ephys_preprocessing/ephys_preprocessing/
```

## Monitoring Jobs

```bash
# View job queue
squeue -u $USER

# View specific job
squeue -j JOB_ID

# View job details
scontrol show job JOB_ID

# Cancel jobs
scancel JOB_ID
```

## Logs

Job logs are written to `logs/`:
- `spikesort_JOBID.out` / `.err` - Single spike sorting job
- `sync_JOBID.out` / `.err` - Single sync job
- `spikesort_ARRAYJOBID_TASKID.out` / `.err` - Array spike sorting
- `sync_ARRAYJOBID_TASKID.out` / `.err` - Array sync

## Troubleshooting

### "Cannot find module" errors
Uncomment and adjust the `module load` lines in the SLURM scripts for your cluster.

### Path not found errors
Check that bind paths in the SLURM scripts match your HPC filesystem layout.

### GPU not detected
- Ensure `--nv` flag is present in `slurm_spikesort.sh`
- Check that your cluster supports GPU access via Singularity
- Verify CUDA modules are loaded

## Notes

- Both services use the **same** Singularity image
- Code is bind-mounted at runtime from `CODE_DIR` — the `.sif` only needs rebuilding when dependencies change
- SLURM handles workflow orchestration (no docker-compose on HPC)
