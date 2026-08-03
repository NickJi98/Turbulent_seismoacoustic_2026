### CM1 Analysis Scripts

This directory contains scripts for analyzing CM1 outputs.

### MATLAB Scripts

The following scripts gather CM1 outputs on cluster.
- cm1_evolution.m: Temporal evolution of diagnostic variables, gathered from each individual diagnostic file.
- cm1_diag_profile.m: Temporal and domain averaged vertical profiles of the final 1-hour of the main run.
- cm1_extract_height.m: Extract variable fields at specific height levels for the dense output run.

The first two scripts can be run with the job submission script `analyze_cm1.sbatch`.
The third one is run with `submit_extract.sh` and `extract_height.sbatch`.

The following scripts visualize results of each individual run.
- init_env.m: Global paths and settings for analyzing outputs of a single CM1 run.
- analyze_evo.m: Visualize temporal evolution of diagnostic variables. This is to check that CM1 run becomes quasi-steady.
- analyze_profiles.m: Visualize vertical profiles of diagnostic variables.
- convec_vel.m: Calculate space-time correlation function in a specified spatial direction.

The following script analyze the parameter-sweeping runs.
- create_db.m: Gather results from the CM1 runs. (e.g. /Data/CM1/param_search/param_search.mat)

### Python Scripts

The following scripts calculate spectral quantities of the turbulent fields.
- calc_spectra_height.py: Pressure and temperature spectra (Welch method).
- calc_spec_vel_height.py: Along-wind and cross-wind spectra (Welch method).
- calc_prms_mpi.py: Calculate RMS pressure amplitude.
- prs_est.py: Obtain mean pressure PSD from the calculated spectra.
- gather_spec.py: Gather spectra from a CM1 run into a single file. (e.g. /Data/CM1/spectra/spectra_*_all.npz)
- plot_spec.py: Plot spectra at different height levels.

The first two scripts can be run with the job submission script `analyze_spec.sbatch`.

### Other Scripts
- gather_diags.sh: Gather all main run outputs into a single directory (e.g. /Data/CM1/param_search/).
                   Then you may run create_db.m to generate a database for plotting.
- submit_extract.sh: Submit a job running cm1_extract_height.m via extract_height.sbatch for each simulation directory.
- submit_py.sh: Submit a job running Python scripts.
- analyze_cm1.sbatch: Submit a job array running cm1_evolution.m and cm1_diag_profile.m for each simulation directory.
- analyze_spec.sbatch: Submit a job array running calc_spec_vel_height.py and calc_spectra_height.py for each directory.
- extract_height.sbatch: Submit a job array running cm1_extract_height.m for different height levels for a single run.

### Cluster Configuration

- Before submitting jobs, modify the SLURM submission scripts to match your cluster settings.
- Adjust partition names, account information, and `ntask-per-node` according to your system.
- SLURM job scripts can be modified to submit individual runs.
