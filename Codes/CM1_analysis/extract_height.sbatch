#!/bin/bash

#SBATCH --job-name=extract_cm1
#SBATCH --output cm1_extract.out.%A
#SBATCH --partition=serc
#SBATCH -c 24
#SBATCH --time=01:00:00
#SBATCH --array=0-4%3

module load matlab

## Header Info
echo;
echo "DATE: $(date)";
start_time=`date +%s`

## Current directory
main_dir=$(pwd)

## Define the list of directories
dir="./V42_Cd22_dense"

## MATLAB script to use
script_dir='/path/to/Codes/CM1_analysis/matlab'
my_script="cm1_extract_height"

## List of height levels
hqs=(15 35 45 55 95)

## Pick hq based on SLURM_ARRAY_TASK_ID
hq=${hqs[$SLURM_ARRAY_TASK_ID]}

echo "Processing hq = $hq in directory: $dir"

## Go into simulation directory
cd "$dir" || { echo "Failed to enter directory: $dir"; exit 1; }

## Copy MATLAB script
cp -u "${script_dir}/${my_script}.m" ./

## Run MATLAB with hq
matlab -batch "hq=${hq}; ${my_script}"

## Return to submission directory
cd "$main_dir"

## Footer Info
echo;
echo "Done...";
end_time=`date +%s`
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES";
echo "Execution time: $(expr $end_time - $start_time) s";
