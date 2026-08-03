#!/bin/bash

#SBATCH --job-name=modeling
#SBATCH --output slurm.out.%j                    
#SBATCH --partition=serc                                        
#SBATCH -c 24
#SBATCH --time=00:40:00

module load matlab/R2023b
module load python/3.9.0 fftw/3.3.10

## CM1 output file
cm1_dir='/your/path/to/cm1/output'
cm1_file='cm1out_prs_15m.nc'
rm -f ./tmp/cm1out_prs.nc
ln -s ${cm1_dir}/${cm1_file} ./tmp/cm1out_prs.nc

## Output directory
out_dir='./outputs/my_run1'
mkdir -p "$out_dir"
escaped_out_dir=$(printf '%s\n' "$out_dir" | sed 's/[&/\]/\\&/g')
sed -i "s|^data_dir = .*|data_dir = '${escaped_out_dir}';|" ./init_env.m
sed -i "s|^les_dir = .*|les_dir = '${escaped_out_dir}'|" ./calc_spectra_mpi.py

## Header Info
echo;
echo "DATE: $(date)";
start_time=`date +%s`

## Run MATLAB script (Change vel_file as needed)
matlab -batch "vel_file='vel_model_fit.csv'; model_main"
matlab -batch model_kernel

## Calculate spectra
python3 calc_spectra_mpi.py

## Footer Info
echo;
echo "Done...";
end_time=`date +%s`
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES";
echo "Execution time: $(expr $end_time - $start_time) s";
