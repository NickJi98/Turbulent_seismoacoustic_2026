#!/bin/bash

#SBATCH --job-name=spectra
#SBATCH --output cm1_spec.out.%A
#SBATCH --partition=serc
#SBATCH --array=0-9
#SBATCH -c 12
#SBATCH --time=00:10:00

module load python/3.9.0 fftw/3.3.10

## Header Info
echo;
echo "DATE: $(date)";
start_time=`date +%s`

## Python script to use
main_dir=$(pwd)
script_dir=$(pwd)
my_scripts=(calc_spec_vel_height.py calc_spectra_height.py plot_spec.py)

## Define the list of directories
directories=($(find *_dense -maxdepth 1 -type d -printf "%f\n" | sort))

## Directory for this job
dir=${directories[$SLURM_ARRAY_TASK_ID]}
echo -e "\nProcessing directory: $dir"

## Run scripts in the chosen directory
cd "$dir" || { echo "Failed to enter directory: $dir"; exit 1; }

for my_script in "${my_scripts[@]}"; do
    echo -e "\nPython script: $my_script"
    cp -u ${script_dir}/${my_script} ./
    python3 ${my_script}
done

## Gather figures
mkdir -p "${main_dir}/figures"
dir_name="${dir%_dense}"
cp "${main_dir}/${dir}/spectra_height/spectra_wind.png" "${main_dir}/figures/spectra_wind_${dir_name}.png"
cp "${main_dir}/${dir}/spectra_height/spectra_tp.png" "${main_dir}/figures/spectra_tp_${dir_name}.png"


## Footer Info
echo;
echo "Done...";
end_time=`date +%s`
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES";
echo "Execution time: $(expr $end_time - $start_time) s";
