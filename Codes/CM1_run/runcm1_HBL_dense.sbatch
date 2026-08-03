#!/bin/bash
#SBATCH --job-name=CM1_HBL                   
#SBATCH --partition=serc
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1          # 1 CPU per task (no multithreading)
#SBATCH --time=12:00:00
#SBATCH --constraint=CPU_MNF:AMD
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
#SBATCH --array=0-24

## Load modules
module load netcdf-fortran/4.5.4

## Directories
cdir=${CM1HOME}
wdir=$(pwd)

## Preparation (Optional)
prep_flag="prep_done.flag"
if [[ $SLURM_ARRAY_TASK_ID -eq 0 && ! -f $prep_flag ]]; then
    echo "Task 0: preparing simulations..."
    bash prepare_sims_dense.sh
    touch $prep_flag
fi

# Wait until preparation is done
while [[ ! -f $prep_flag ]]; do
    sleep 5
done

## Directories for this job
# Example: sim_dirs=(dir1_dense dir2_dense dir3_dense)
sim_dirs=( $(find . -maxdepth 1 -type d -name "V*_Cd*_dense" -printf "%f\n" | sort) )
simdir=${sim_dirs[$SLURM_ARRAY_TASK_ID]}
echo -e "\nTask $SLURM_ARRAY_TASK_ID processing directory: $simdir"

#### Go to ${wdir} ####
cd ${wdir}

## Copy executable and data files
cp $cdir/run/LANDUSE.TBL     	$simdir/
cp $cdir/run/cm1.exe    	$simdir/
cp $cdir/run/onefile.F      	$simdir/
cp $cdir/run/RRTMG_LW_DATA 	$simdir/
cp $cdir/run/RRTMG_SW_DATA   	$simdir/

#### Go to ${simdir} ####
cd $simdir

## Header Info
echo "DATE: $(date)"
start_time=`date +%s`

## Run simulation
mpirun cm1.exe > OUT &
wait

## Footer Info
echo
echo "Done..."
end_time=`date +%s` 
echo "DATE: $(date), CPUS_PER_TASK: $SLURM_CPUS_PER_TASK, NNODES: $SLURM_NNODES"

## Move output file
mv ${wdir}/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out \
   ${wdir}/${simdir}/runtime_${SLURM_ARRAY_TASK_ID}.info

## Check error output
error_file="${wdir}/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

# Check if the error output file is empty
if [ ! -s "$error_file" ]; then
    echo "The error output file is empty. Deleting..."
    rm "$error_file"
else
    echo "The error output file is not empty."
fi
