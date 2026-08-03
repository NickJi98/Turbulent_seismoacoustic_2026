### Quasi-static elastic modeling

This directory contains scripts for quasi-static elastic modeling of vertical displacement under turbulent pressure loading.

For using the propagator matrix method for quasi-static modeling, you may also check the following GitHub repository for details: https://github.com/NickJi98/Quasi-static-deformation

### Scripts

You can use the SLURM job submission scripts `submit_model.sh` to perform the modeling (or `submit_model_array.sh` for modeling with several velocity profiles). Set the location of CM1 pressure field and velocity model file before running. The outputs will be placed under /outputs/.

The velocity model file `vel_model_fit.csv' is used for Figure 3, while `vel_model_site.csv' is from the survey site for comparison. They are copied from /Data/vel_model.

`model_main.m` performs the quasi-static modeling. `model_kernel.m` performs modeling in perturbed velocity models. Each layer is perturbed by 1% of shear-wave velocity (Vs), with Vp and rho modified accordingly.

`calc_spectra_mpi.py` calculates the spectra of all the modeling results.
