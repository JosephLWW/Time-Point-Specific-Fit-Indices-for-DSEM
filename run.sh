#!/bin/bash
#SBATCH --job-name=hdi-pipeline
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1

#SBATCH --array=75-124
#SBATCH --job-name=rslurm_multi_array
#SBATCH --output=slurm_%a.out
#SBATCH --time=24:00:00
#SBATCH --partition=single
#SBATCH --mem-per-cpu=500

module purge
module load math/R/4.1.2




# Define Variables
SCRIPT_DIR="scripts"
LOG_DIR="logs"
DATA_DIR="data/raw"

# Start Message
echo "=========================================="
echo "Starting HDI Pipeline at: $(date)"
echo "=========================================="

# Check/Create Directories
if [ ! -d "$LOG_DIR" ]; then
    echo "Creating logs directory..."
    mkdir -p "$LOG_DIR"
fi


# Execute
echo "Running..."
Rscript main.R > "$LOG_DIR/log_$SLURM_ARRAY_TASK_ID.txt" 2>&1

# 5. Check for Success (Exit Code Analysis)
EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo " completed successfully."
else
    echo "ERROR: exit code $EXIT_CODE."
    exit 1  # Stop the pipeline here
fi

# 7. Final Message
echo "=========================================="
echo "Pipeline finished successfully at: $(date)."
echo "=========================================="