#!/usr/bin/env bash
#SBATCH --job-name=TaskFarm
#SBATCH --partition=modi_HPPC
#SBATCH --exclusive
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

echo "Host: $(hostname)"
echo "NP (SLURM_NTASKS): $SLURM_NTASKS"
echo "Workers: $((SLURM_NTASKS-1))"

IMAGE="$HOME/modi_images/hpc-notebook-latest.sif"

# Run inside container (Slurm launches MPI ranks; no mpiexec needed)
srun --cpu-bind=cores singularity exec "$IMAGE" ./task_farm_HEP
