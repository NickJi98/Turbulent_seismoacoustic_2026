#!/usr/bin/env python3

"""
Spectral analysis of CM1 modeling results
(Input horizontal wind fields at specific height levels)

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

# Extract height levels
# height_list = [5,15,25,35,45,55,75,95,145,195]
height_list = [15]

# Welch parameters
dt = 1.0
nperseg = 1024
noverlap = int(nperseg * 0.5)

for height in height_list:

    # Read variable field
    _u = nc.Dataset(les_dir + f'cm1out_u_{height}m.nc', 'r')
    u_field = _u.variables['u'][:].transpose(2, 1, 0)
    _v = nc.Dataset(les_dir + f'cm1out_v_{height}m.nc', 'r')
    v_field = _v.variables['v'][:].transpose(2, 1, 0)
    Nx, Ny = v_field.shape[0], u_field.shape[1]

    # Frequency samples
    _, freqs = fspec.my_welch(np.squeeze(u_field[1,1,:]), dt, nperseg)
    Nf = len(freqs)

    # Interpolate to common grid (256×256)
    u_ = 0.5 * (u_field[:-1, :, :] + u_field[1:, :, :])
    v_ = 0.5 * (v_field[:, :-1, :] + v_field[:, 1:, :])

    # Mean wind direction
    u_avg, v_avg = np.mean(u_), np.mean(v_)
    theta = np.arctan2(v_avg, u_avg)
    print(f'Height {height} m. Mean wind theta {np.rad2deg(theta)} deg')

    # Rotation
    u_along = u_ * np.cos(theta) + v_ * np.sin(theta)
    u_cross = -u_ * np.sin(theta) + v_ * np.cos(theta)

    # Perturbation
    u_along = u_along - np.mean(u_along)
    u_cross = u_cross - np.mean(u_cross)

    # Parallel implementation
    u_along_shm = shared_memory.SharedMemory(create=True, size=u_along.nbytes)
    u_along_shared = np.ndarray(u_along.shape, dtype=u_along.dtype, buffer=u_along_shm.buf)
    u_along_shared[:] = u_along[:]
    u_cross_shm = shared_memory.SharedMemory(create=True, size=u_cross.nbytes)
    u_cross_shared = np.ndarray(u_cross.shape, dtype=u_cross.dtype, buffer=u_cross_shm.buf)
    u_cross_shared[:] = u_cross[:]

    # Compute PSDs
    def _func(i, j, u_shape, dt, nperseg, noverlap):

        # Access shared memory
        u1 = np.ndarray(u_shape, dtype=np.float32, buffer=u_along_shm.buf)
        u2 = np.ndarray(u_shape, dtype=np.float32, buffer=u_cross_shm.buf)

        # Compute PSDs and transfer function
        psd1_, _ = fspec.my_welch(np.squeeze(u1[i,j,:]), dt, nperseg, noverlap, ave_method='mean')
        psd2_, _ = fspec.my_welch(np.squeeze(u2[i,j,:]), dt, nperseg, noverlap, ave_method='mean')
        return psd1_, psd2_

    # Parallelization
    with Pool(processes=n_cpus) as pool:
        results = pool.starmap(_func, [(i, j, u_along.shape, dt, nperseg, noverlap)
                                       for i in range(Nx) for j in range(Ny)])
    # Reshape result into array
    psd1_les, psd2_les = np.zeros((Nx,Ny,Nf)), np.zeros((Nx,Ny,Nf))
    for idx, (i, j) in enumerate([ (i, j) for i in range(Nx) for j in range(Ny) ]):
        psd1_les[i, j, :], psd2_les[i, j, :] = results[idx]

    # Clean up shared memory
    u_along_shm.close()
    u_along_shm.unlink()
    u_cross_shm.close()
    u_cross_shm.unlink()

    # Median and interquartile range
    psd_u1, psd_u2 = np.zeros((3,Nf)), np.zeros((3,Nf))
    psd_u1[0, :] = np.median(psd1_les, axis=(0,1))
    psd_u1[1, :] = np.percentile(psd1_les, 25, axis=(0,1))
    psd_u1[2, :] = np.percentile(psd1_les, 75, axis=(0,1))
    psd_u2[0, :] = np.median(psd2_les, axis=(0,1))
    psd_u2[1, :] = np.percentile(psd2_les, 25, axis=(0,1))
    psd_u2[2, :] = np.percentile(psd2_les, 75, axis=(0,1))

    # Save median spectra
    np.savez(os.path.join(spec_dir, f'spectra_wind_{height:d}m.npz'), 
             freqs=freqs, psd_u1=psd_u1, psd_u2=psd_u2)
