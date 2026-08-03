#!/bin/bash

# Slurm job submission script
sbatch_file="extract_height.sbatch"

# Loop through all simulation directories
dirs=($(find V*_Cd*_dense -maxdepth 1 -type d -printf "%f\n" | sort))

for dir in "${dirs[@]}"; do
    if [ -d "$dir" ]; then
        echo "Submitting job for simulation directory: ${dir}"

        # Create a unique job script for each directory
        job_file="${dir}_extract.sbatch"
        sed "s|^dir=.*|dir='${dir}'|" "$sbatch_file" > "$job_file"

        # Submit the job
        sbatch "$job_file"
    fi
done
