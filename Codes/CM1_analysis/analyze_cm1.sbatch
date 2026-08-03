#!/bin/bash

#SBATCH --job-name=analyze_cm1
#SBATCH --output cm1_analysis_%a.out                    
#SBATCH --partition=serc                                        
#SBATCH -c 12
#SBATCH --time=00:30:00
#SBATCH --array=0-24

module load matlab

## Header Info
echo;
echo "DATE: $(date)";
start_time=`date +%s`

## Current directory
main_dir=$(pwd)

## Define the list of directories
directories=( $(find "${main_dir}" -maxdepth 1 -type d -name "V*_Cd*" ! -name "*_dense" -printf "%f\n" | sort) )

## Directory for this job
dir=${directories[$SLURM_ARRAY_TASK_ID]}
echo -e "\nProcessing directory: $dir"

## MATLAB script to use
script_dir='/path/to/my_matlab_scripts'
my_scripts=(cm1_evolution cm1_diag_profile)

## Loop over the directories
for my_script in "${my_scripts[@]}"; do
    echo -e "\n"
    echo "MATLAB script: $my_script"
    
    ## Enter the directory
    cd "$dir" || { echo "Failed to enter directory: $dir"; exit 1; }
    
    ## Copy MATLAB script
    cp -u ${script_dir}/${my_script}.m ./
    
    ## Run MATLAB script
    matlab -batch ${my_script}
    
    ## Return to submission directory
    cd $main_dir
done

## Copy files
cd "$dir" || { echo "Failed to enter directory: $dir"; exit 1; }
cp namelist.input ./Data/
cp cm1out_stats.nc ./Data/
cp cm1out_diag_000001.nc ./Data/base_state.nc
cd $main_dir

## Footer Info
echo;
echo "Done...";
end_time=`date +%s`
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES";
echo "Execution time: $(expr $end_time - $start_time) s";
