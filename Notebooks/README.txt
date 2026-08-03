### Jupyter Notebooks

This directory contains the Jupyter Notebooks performing analyses and reproducing figures of the paper.

Before running the notebooks:
- Build the python environment from the top-level `environment.yml` (see top-level README.txt).
- Download the Datasets (https://doi.org/10.25740/zk709vr3334) and place the `Data/` folder at the top level of this archive, next to `Notebooks/`.

Individual figures are saved to /Figures/Panels/ and /Figures_Sup/Panels/ (the figure saving lines are commented out by default).

### Common Setup

- _common.py: Common paths, python packages and helper functions. Loaded by every notebook via `from _common import *`.

This file derives all data paths from its own location (no manual path editing needed if `Data/` is at the top level). It also appends `Codes/python_packages` to the python path and applies the shared plotting style `Codes/python_packages/matplotlibrc`.

### Main Notebooks

- Maps.ipynb: Maps plotted on top of the satellite image.
- Spectrogram.ipynb: Waveforms and spectrograms of pressure and seismic data.
- Wind_Analysis.ipynb: Wind speeds from ASOS, ERA5 and NEXRAD VAD profiles, and boundary layer height scales.
- CM1_inputs.ipynb: Hurricane datasets (radar, dropsondes, reanalysis) and vertical profiles for the CM1 setup.
- Modeling.ipynb: Visualizing LES pressure fields and elastic modeling results.
- Transfer.ipynb: Transfer function between pressure and ground displacement.
- Turb_Analysis.ipynb: Turbulent dissipation rates estimated from wind/pressure spectra and structure functions.

### Supplementary Notebooks under calc_results/

These notebooks download seismometer and infrasound data and perform the heavier calculations to generate intermediate results. The outputs are already included in the Datasets, so they only need to be re-run to reproduce those files from raw data.

- calc_obs_spec.ipynb: Observed wind-tower and seismic/infrasound spectra (outputs under /Data/Isaac/obs_spectra.npz).
- calc_spectrogram.ipynb: Wavelet (fCWT) spectrograms, cross-spectra and coherence (outputs under /Data/TA/spectrogram/).
- calc_turb.ipynb: Turbulent pressure PSD, structure functions and RMS amplitudes (outputs under /Data/TA/turbulence/).

### Notes

- The fCWT package is only required for the wavelet-based notebooks (calc_spectrogram.ipynb, calc_obs_spec.ipynb). The other notebooks run without it (a warning is printed if fCWT is missing).
- The raw tower data are from the Florida Coastal Monitoring Program and are available upon request (see Acknowledgments of the paper).
