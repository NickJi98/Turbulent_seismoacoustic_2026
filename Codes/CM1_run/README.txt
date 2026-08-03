### CM1 Simulation Setup Scripts

This directory contains scripts for preparing and running parameter-sweep simulations with CM1.

### Files

- prepare_sims.sh: Creates simulation directories for different input gradient wind speed (V) and drag coefficient (Cd) combinations.
- prepare_sims_dense.sh: After the main runs, creates directories for 1 Hz output simulations.

- SLURM job scripts (`runcm1_HBL.sbatch`, `runcm1_HBL_dense.sbatch`): Submit job arrays for each prepared directory.

### Input Files

- The files `input_sounding` and `namelist.input` (`namelist.input_dense`) are the same as those in `/Data/CM1/inputs`.
- Keep them synchronized if the base configuration changes.

### Main Run vs. Dense Output Run

- Main runs ensure that the HBL simulations become quasi-steady. It runs 6-hour simulation.
- Dense output run writes output every 1 s to obtain the turbulent fields for spectral analysis. It restarts from the final snapshot of the main run, and simulates for an additional 1 hour.

### Cluster Configuration

- Before submitting jobs, modify the SLURM submission scripts to match your cluster settings.
- Adjust partition names, account information, and `ntask-per-node` according to your system.
- Ensure that the `ppnode` parameter (MPI processes per node) in `namelist.input` is consistent with your cluster settings.
- SLURM job scripts can be modified to submit individual runs.
