#!/bin/bash
#SBATCH --job-name=dsem_sim
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --array=1-324%64
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

# Clean modules and load R
module purge
# To check R version: module avail R
module load math/R

if [ -z "$SLURM_ARRAY_JOB_ID" ]; then
    export RUN_ID="${SLURM_JOB_ID}"
else
    export RUN_ID="${SLURM_ARRAY_JOB_ID}"
fi

echo "Master RUN_ID is: $RUN_ID"

LOG_DIR="logs"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi

echo "=========================================="
echo "Starting DSEM Simulation"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "=========================================="

# Ejecutar R (Sin redefinir RUN_ID)
Rscript main.R

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS"
    rm "logs/sim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
    rm "logs/sim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
else
    echo "ERROR"
    exit 1
fi