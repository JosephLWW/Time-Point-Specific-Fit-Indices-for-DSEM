#!/bin/bash
#SBATCH --job-name=dsem_sim
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --array=1-324%64
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=256000
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

# Clean modules and load R
module purge
# To check R version: module avail R
module load math/R

# Set a global Run ID for this execution
export RUN_ID="${SLURM_ARRAY_JOB_ID}_$(date +%Y%m%d_%H%M%S)"
echo "Master RUN_ID is: $RUN_ID"

# 2. Definir Directorio log
LOG_DIR="logs"

# 3. Crear directorio de logs si no existe
if [ ! -d "$LOG_DIR" ]; then
    echo "Creating logs directory..."
    mkdir -p "$LOG_DIR"
fi


# 4. Mensaje de Inicio
echo "=========================================="
echo "Starting DSEM Simulation at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $HOSTNAME"
echo "=========================================="

# 5. Ejecutar R
RUN_ID="${SLURM_JOB_ID}_$(date +%Y%m%d_%H%M%S)"
export RUN_ID
Rscript main.R

# 6. Verificación de Salida
EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: R script finished correctly."
else
    echo "ERROR: R script failed with exit code $EXIT_CODE."
    exit 1
fi

echo "=========================================="
echo "Finished at: $(date)"
echo "=========================================="
