# Turbulent seismoacoustic imprints during a hurricane landfall

Codes for the following study:

> Ji, Q., Dey, I., and Dunham, E. (2026). Turbulent seismoacoustic imprints during a hurricane landfall.

The Datasets are available at <https://doi.org/10.25740/zk709vr3334>

Download the Datasets and place the `Data/` folder at the top level of this directory (next to `Codes/` and `Notebooks/`). The Jupyter Notebooks assume this location; for the MATLAB scripts, set your data directories in the respective `init_env.m`. See [Data/README.txt](Data/README.txt) for the dataset contents.

## Directory Structure

- `Codes/python_packages/`: Python packages (seismic data acquisition, Fourier spectral analysis, atmospheric datasets) and the `matplotlibrc` plotting style.
- `Codes/CM1_run/`: Preparing and running CM1 large-eddy simulations ([README](Codes/CM1_run/README.txt)).
- `Codes/CM1_analysis/`: Analyzing CM1 outputs with MATLAB and python ([README](Codes/CM1_analysis/README.txt)).
- `Codes/elastic_modeling/`: Quasi-static elastic modeling of vertical displacement ([README](Codes/elastic_modeling/README.txt)).
- `Notebooks/`: Jupyter Notebooks performing analyses and reproducing the figures ([README](Notebooks/README.txt)).
- `Figures/`, `Figures_Sup/`: Main and supplementary figures, with individual figure panels under `Panels/`.
- `fcwt_install/`: Patched source files and installation script for the fCWT package (might fail in Mac environment).

## Python Environment

Build the conda environment, then install fCWT (next section):

```bash
conda env create -f environment.yml
conda activate seismo-isaac
bash fcwt_install/install_fcwt.sh   # optional
```

The environment name is set by the `name:` field in `environment.yml` (default: `seismo-isaac`); feel free to change it.

- Tested python versions: 3.11 (Mac, Apple M3 chip), 3.9 / 3.12 (Linux HPC cluster).
- To run the python scripts under `Codes/CM1_analysis/python/` and `Codes/elastic_modeling/` on a cluster, add `Codes/python_packages` to your PYTHONPATH.

## Install fCWT Package

> Arts, L.P.A., van den Broek, E.L. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time-frequency analysis. Nat Comput Sci 2, 47-58 (2022).

Package source: <https://github.com/fastlib/fCWT>

fCWT is built from source because PyPI provides no prebuilt wheels for Apple Silicon Mac, Linux or Windows. The script `fcwt_install/install_fcwt.sh` automates the steps below. fCWT is only required for the wavelet-based notebooks (see [Notebooks/README.txt](Notebooks/README.txt)); the other notebooks run without it.

Last tested: 2025-05-06

Requirement: python >= 3.7, numpy >= 1.14.5, matplotlib

### Mac Environment (Apple M3 chip, python=3.11)

1. Package source

   ```bash
   git clone https://github.com/fastlib/fCWT.git
   ```

2. Replace / Add several files

   Replace / add the source files with those provided under `fcwt_install/`. The modifications to each original file are explained.

   - `fCWT/setup.py`: modify `comp_args`, `library_dirs`, `include_dirs`
   - `fCWT/src/fcwt/boilerplate.py`: fix a bug in the function `plot`
   - `fCWT/src/fcwt/fcwt.cpp`: add `#include <cassert>`
   - `fCWT/libs/libomp.dylib`: add dynamic library from `libomp`

3. Install and test

   ```bash
   pip install .
   ```

   Then you can try the python example from GitHub.

### Linux Environment (e.g. HPC cluster, python=3.9/3.12)

For HPC cluster, after loading python, you may also need to load fftw (e.g. 3.3.10).

Replace `boilerplate.py`, `fcwt.cpp`. For `fcwt.cpp` you may also need to add `#include <cstring>`.

## Jupyter Notebooks

The notebooks under `Notebooks/` reproduce the observational analysis and modeling figures. Common paths, packages and helper functions are loaded from `Notebooks/_common.py`, which also applies the plotting style `Codes/python_packages/matplotlibrc`.

See [Notebooks/README.txt](Notebooks/README.txt) for the list of notebooks and details.

## Python Packages (Codes/python_packages)

- `seismic/`: Seismic station data acquisition and plotting (obspy-based).
- `Fourier/`: Fourier and wavelet spectral analysis (Welch, CWT/fCWT, MISO coherence analysis, PSD-PDF).
- `TC_analysis/`: Analysis of atmospheric / hurricane datasets and CM1 outputs.
- `matplotlibrc`: Plotting style, applied by `Notebooks/_common.py`.

## MATLAB Scripts

The Parallel Computing Toolbox is used for elastic modeling.

- Quasi-static elastic response modeling: see [Codes/elastic_modeling/README.txt](Codes/elastic_modeling/README.txt)
- Analysis of LES results from CM1: see [Codes/CM1_analysis/README.txt](Codes/CM1_analysis/README.txt)

Make sure to update your data directories in `init_env.m` before running the scripts. The raw outputs from CM1 are already processed and included in the Datasets (`Data/CM1/outputs/`).

## CM1 Simulations

- Setup and job submission for the CM1 runs: see [Codes/CM1_run/README.txt](Codes/CM1_run/README.txt)
- Before submitting jobs, modify the SLURM scripts (under `Codes/CM1_run/`, `Codes/CM1_analysis/` and `Codes/elastic_modeling/`) to match your cluster settings.
