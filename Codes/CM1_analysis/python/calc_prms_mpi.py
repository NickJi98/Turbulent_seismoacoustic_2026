#!/usr/bin/env python3

"""
Spectral analysis of CM1 modeling results

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np
import netCDF4 as nc
import pandas as pd
from obspy.signal import filter
from multiprocessing import shared_memory, Pool


# Number of CPUs
n_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"Using {n_cpus} CPUs")

# Directory for CM1 modeling outputs
les_dir = '/path/to/simulations'    # MODIFY THIS PATH

# Get all simulation directories
sim_dirs = sorted([d for d in os.listdir(les_dir) 
			    if os.path.isdir(os.path.join(les_dir, d)) and d.endswith("dense")])
nc_file = 'cm1out_prs_15m.nc'

# Store RMS pressure
p_rms_arr = []

for sim_dir in sim_dirs:
    print(f"Processing directory: {sim_dir}")

    # Pressure field
    les_result = nc.Dataset(os.path.join(les_dir, sim_dir, nc_file), 'r')
    psfc = les_result.variables['prs'][:].transpose(2, 1, 0)
    _shape = psfc.shape
    Nx, Ny, Nt = _shape[0], _shape[1], _shape[2]

    # Perturbation pressure
    p_mean = np.mean(psfc)
    pp = psfc - p_mean
    pp = pp - np.mean(pp, axis=(0,1))

    ### Parallel Implementation ###

    # Create shared memory
    pp_shm = shared_memory.SharedMemory(create=True, size=pp.nbytes)

    # Copy data into shared memory
    pp_shared = np.ndarray(pp.shape, dtype=pp.dtype, buffer=pp_shm.buf)
    pp_shared[:] = pp[:]

    # Filter
    def _func(i, j, pp_shape):

        # Access shared memory
        pp = np.ndarray(pp_shape, dtype=np.float32, buffer=pp_shm.buf)

        # Filtering
        pp_filt = filter.bandpass(np.squeeze(pp[i, j, :]), freqmin=0.005, freqmax=0.2, df=1)
        
        return pp_filt

    # Parallelization
    with Pool(processes=n_cpus) as pool:
        results = pool.starmap(_func, [(i, j, pp.shape) for i in range(Nx) for j in range(Ny)])

    # Reshape results into arrays
    pp_filt = np.zeros((Nx,Ny,Nt))
    for idx, (i, j) in enumerate([(i, j) for i in range(Nx) for j in range(Ny)]):
        pp_filt[i, j, :] = results[idx]

    # Clean up shared memory
    pp_shm.close()
    pp_shm.unlink()

    # RMS pressure
    p_rms = np.std(pp_filt)
    print(f'RMS pressure: {p_rms:.2f} Pa')
    p_rms_arr.append(p_rms)

# Save to CSV file
df = pd.DataFrame({"case": sim_dirs, "p_rms": p_rms_arr})
df.to_csv("output_5-200s.csv", index=False)
