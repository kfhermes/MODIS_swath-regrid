#!/bin/bash

#SBATCH --job-name=regrid_modis_aqua   # Job name
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=swift
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem=6000mb                # Job Memory
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.log    # Standard output and error log
#SBATCH --array=2000-2024           # Array range
pwd; hostname; date

#!/bin/bash
module load jaspy

YEAR=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")

# Define the base directory
BASE_DIR="/neodc/modis/data/MOD04_L2/collection61"
ODIR="directory/where/you/want/to/save/the/output"

# Loop through months
for MONTH in {01..12}; do
    # Loop through days
    for DAY in {01..31}; do
        # Define the directory path
        IDIR="${BASE_DIR}/${YEAR}/${MONTH}/${DAY}"
        # Check if the directory exists
        if [ -d "$IDIR" ]; then
            echo "Processing $IDIR"
            python compute_dod_regrid_xesmf.py "${IDIR}" "${ODIR}"
        fi
    done
done