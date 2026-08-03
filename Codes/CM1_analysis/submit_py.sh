#!/bin/bash

#SBATCH --job-name=spectra
#SBATCH --output slurm.out.%j                    
#SBATCH --partition=serc                                        
#SBATCH -c 12
##SBATCH --mem-per-cpu=2GB
#SBATCH --time=00:60:00
##SBATCH --dependency=afterok:1527862

module load python/3.9.0 fftw/3.3.10

## Header Info
echo;
echo "DATE: $(date)";
start_time=`date +%s`

## Run Python script
python3 calc_prms_mpi.py

## Footer Info
echo;
echo "Done...";
end_time=`date +%s`
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES";
echo "Execution time: $(expr $end_time - $start_time) s";