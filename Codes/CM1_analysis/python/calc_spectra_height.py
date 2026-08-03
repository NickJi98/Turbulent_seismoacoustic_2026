#!/usr/bin/env python3

"""
Spectral analysis of CM1 modeling results
(Input variable fields at specific height levels)

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np
import netCDF4 as nc
from multiprocessing import shared_memory, Pool

import Fourier.spectrum as fspec


# Number of CPUs
n_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"Using {n_cpus} CPUs")

# Directory for CM1 results
les_dir = './'
spec_dir = os.path.join(les_dir, 'spectra_height')
os.makedirs(spec_dir, exist_ok=True)

# Variable names
var_list = ('prs', 'th')

# Extract height levels
# height_list = [5,15,25,35,45,55,75,95,145,195]
height_list = [15]

# Welch parameters
dt = 1.0
nperseg = 1024
noverlap = int(nperseg * 0.5)

for varname in var_list:
    for height in height_list:

        # Read variable field
        nc_file = f'cm1out_{varname}_{height:d}m.nc'
        data_ = nc.Dataset(les_dir + nc_file, 'r')
        field = data_.variables[varname][:].transpose(2, 1, 0)
        Nx, Ny = field.shape[0], field.shape[1]

        # Perturbation
        field = field - np.mean(field)

        # Frequency samples
        _, freqs = fspec.my_welch(np.squeeze(field[1,1,:]), dt, nperseg)
        Nf = len(freqs)

        # Parallel implementation
        field_shm = shared_memory.SharedMemory(create=True, size=field.nbytes)
        field_shared = np.ndarray(field.shape, dtype=field.dtype, buffer=field_shm.buf)
        field_shared[:] = field[:]

        # Compute PSDs
        def _func(i, j, field_shape, dt, nperseg, noverlap):

            # Access shared memory
            field = np.ndarray(field_shape, dtype=np.float32, buffer=field_shm.buf)

            # Compute PSDs and transfer function
            psd_, _ = fspec.my_welch(np.squeeze(field[i,j,:]), dt, nperseg, noverlap, ave_method='mean')  
            return psd_

        # Parallelization
        with Pool(processes=n_cpus) as pool:
            results = pool.starmap(_func, [(i, j, field.shape, dt, nperseg, noverlap)
                                           for i in range(Nx) for j in range(Ny)])
        # Reshape result into array
        psd_les = np.zeros((Nx,Ny,Nf))
        for idx, (i, j) in enumerate([(i, j) for i in range(Nx) for j in range(Ny)]):
            psd_les[i, j, :] = results[idx]

        # Clean up shared memory
        field_shm.close()
        field_shm.unlink()

        # Median and interquartile range
        psd = np.zeros((3,Nf))
        psd[0, :] = np.median(psd_les, axis=(0,1))
        psd[1, :] = np.percentile(psd_les, 25, axis=(0,1))
        psd[2, :] = np.percentile(psd_les, 75, axis=(0,1))

        # Save median spectra
        np.savez(os.path.join(spec_dir, f'spectra_{varname}_{height:d}m.npz'), 
                 freqs=freqs, psd=psd)
