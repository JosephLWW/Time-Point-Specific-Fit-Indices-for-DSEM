#!/bin/bash
#SBATCH --job-name=dsem_sim
#SBATCH --partition=cpu              # <--- CORREGIDO (antes era 'single')
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # 1 CPU por simulación es suficiente usualmente
#SBATCH --mem-per-cpu=4G             # <--- AUMENTADO (500MB es muy poco para Stan/DSEM)
#SBATCH --time=04:00:00              # Dale tiempo suficiente a Stan (ajusta según necesidad)
#SBATCH --array=75-124               # Tu rango de arrays
#SBATCH --output=logs/sim_%A_%a.out  # Logs organizados (%A=jobID, %a=arrayID)
#SBATCH --error=logs/sim_%A_%a.err

# 1. Limpiar y Cargar Módulos
module purge
# Asegúrate que esta versión de R existe con: module avail R
module load math/R/4.1.2 

# 2. Definir Directorios (Basado en la estructura que creamos)
# Usamos $SLURM_SUBMIT_DIR para asegurar que estamos en la raíz del proyecto
PROJECT_ROOT=$SLURM_SUBMIT_DIR
SCRIPT_PATH="$PROJECT_ROOT/scripts/01_model_fitting/main.R"
LOG_DIR="$PROJECT_ROOT/logs"

# 3. Crear directorio de logs si no existe
mkdir -p "$LOG_DIR"

# 4. Mensaje de Inicio
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $HOSTNAME"
echo "Starting DSEM Simulation at: $(date)"
echo "Running script: $SCRIPT_PATH"
echo "=========================================="

# 5. Ejecutar R
# Pasamos el ID del array como argumento al script de R
Rscript "$SCRIPT_PATH" "$SLURM_ARRAY_TASK_ID"

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
