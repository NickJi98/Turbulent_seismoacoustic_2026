#!/usr/bin/env python3

"""
Spectral analysis of CM1 modeling results

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np
import netCDF4 as nc
from multiprocessing import shared_memory, Pool

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

import Fourier.spectrum as fspec
from scipy.stats import circmean, circstd


# Number of CPUs
n_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
print(f"Using {n_cpus} CPUs")

# Directory for CM1 modeling outputs
les_dir = './outputs/your_directory/'

# Loop over modeling results
nc_files = ['LES_ref.nc', 'LES_kernel_1.nc', 'LES_kernel_2.nc',
            'LES_kernel_3.nc', 'LES_kernel_4.nc', 
            'LES_kernel_5.nc', 'LES_kernel_6.nc', 
            'LES_kernel_7.nc']
            
for nc_file in nc_files:
    les_result = nc.Dataset(os.path.join(les_dir, nc_file), 'r')
    pp = les_result.variables['pp'][:].transpose(2, 1, 0)
    uz = les_result.variables['uz_layer'][:].transpose(2, 1, 0)
    
    # Welch parameters
    dt = 1.0
    nperseg = 1024
    noverlap = int(nperseg * 0.5)
    Nx = pp.shape[0]
    
    # Frequency samples
    _, freqs_les = fspec.my_welch(np.squeeze(pp[1,1,:]), dt, nperseg)
    Nf = len(freqs_les)
    
    ### Parallel Implementation ###
    
    # Create shared memory
    pp_shm = shared_memory.SharedMemory(create=True, size=pp.nbytes)
    uz_shm = shared_memory.SharedMemory(create=True, size=uz.nbytes)
    
    # Copy data into shared memory
    pp_shared = np.ndarray(pp.shape, dtype=pp.dtype, buffer=pp_shm.buf)
    uz_shared = np.ndarray(uz.shape, dtype=uz.dtype, buffer=uz_shm.buf)
    pp_shared[:] = pp[:]
    uz_shared[:] = uz[:]
    
    # Compute PSDs and transfer function
    def _func(i, j, pp_shape, uz_shape, dt, nperseg, noverlap):
    
        # Access shared memory
        pp = np.ndarray(pp_shape, dtype=np.float32, buffer=pp_shm.buf)
        uz = np.ndarray(uz_shape, dtype=np.float32, buffer=uz_shm.buf)
        
        # Compute PSDs and transfer function
        psd_P, _ = fspec.my_welch(np.squeeze(pp[i, j, :]), dt, nperseg, noverlap, ave_method='mean')
        psd_Z, _ = fspec.my_welch(np.squeeze(uz[i, j, :]), dt, nperseg, noverlap, ave_method='mean')
        csd_PZ, _ = fspec.my_welch_csd(np.squeeze(pp[i, j, :]), np.squeeze(uz[i, j, :]), 
                                                  dt, nperseg, noverlap, ave_method='mean')
        tf = csd_PZ / psd_P
        
        return psd_P, psd_Z, tf
    
    # Parallelization
    with Pool(processes=n_cpus) as pool:
        results = pool.starmap(_func, [(i, j, pp.shape, uz.shape, dt, nperseg, noverlap)
                                       for i in range(Nx) for j in range(Nx)])
    
    # Reshape results into arrays
    psd_P_les, psd_Z_les, tf_les = np.zeros((Nx,Nx,Nf)), np.zeros((Nx,Nx,Nf)), np.zeros((Nx,Nx,Nf))
    for idx, (i, j) in enumerate([(i, j) for i in range(Nx) for j in range(Nx)]):
        psd_P_les[i, j, :], psd_Z_les[i, j, :], tf_les[i, j, :] = results[idx]
    
    # Clean up shared memory
    pp_shm.close()
    uz_shm.close()
    pp_shm.unlink()
    uz_shm.unlink()
    
    # Median and interquartile range
    psd_P, psd_Z, tf_PZ, ph_PZ = np.zeros((3,Nf)), np.zeros((3,Nf)), np.zeros((3,Nf)), np.zeros((2,Nf))
    psd_P[0, :] = np.median(psd_P_les, axis=(0,1))
    psd_P[1, :] = np.percentile(psd_P_les, 25, axis=(0,1))
    psd_P[2, :] = np.percentile(psd_P_les, 75, axis=(0,1))
    psd_Z[0, :] = np.median(psd_Z_les, axis=(0,1))
    psd_Z[1, :] = np.percentile(psd_Z_les, 25, axis=(0,1))
    psd_Z[2, :] = np.percentile(psd_Z_les, 75, axis=(0,1))
    tf_PZ[0, :] = np.median(np.abs(tf_les), axis=(0,1))
    tf_PZ[1, :] = np.percentile(np.abs(tf_les), 25, axis=(0,1))
    tf_PZ[2, :] = np.percentile(np.abs(tf_les), 75, axis=(0,1))
    ph_PZ[0, :] = circmean(np.angle(tf_les), axis=(0,1))
    ph_PZ[1, :] = circstd(np.angle(tf_les), axis=(0,1))
    
    # Save median spectra: Reference model
    basename = os.path.splitext(nc_file)[0].replace('LES_', '')
    np.savez(os.path.join(les_dir, f'spectra_{basename}.npz'), 
             freqs=freqs_les, psd_P=psd_P, psd_Z=psd_Z, tf_PZ=tf_PZ, ph_PZ=ph_PZ)

# Plot pressure spectra
ref_data = np.load(os.path.join(les_dir, f'spectra_ref.npz'))
freqs, psd_P = ref_data['freqs'], ref_data['psd_P']

fig, axes = plt.subplots(1, 2, figsize=(5, 2.5), sharex=True)
axes[0].loglog(1/freqs, psd_P[0], 'r-', lw=2)
axes[0].fill_between(1/freqs, psd_P[1], psd_P[2], color='red', alpha=0.2)

# Vertical displacement
axes[1].loglog(1/freqs, psd_Z[0], 'r-', label='Modeling Result', lw=2)
axes[1].fill_between(1/freqs, psd_Z[1], psd_Z[2], color='red', alpha=0.2)

# Modification
axes[0].set_xlabel('Period (s)')
axes[0].set_ylabel('PSD (Pa$^2$/Hz)')
axes[0].set_xlim(1/np.flip([5e-3, 0.25]))
axes[0].set_ylim([1e1, 2e3])
axes[0].grid(which='both')
axes[0].set_title('Surface Pressure')
axes[1].set_xlabel('Period (s)')
axes[1].set_ylabel('PSD ($\mu$m$^2$/Hz)')
axes[1].set_ylim([1e0, 5e2])
axes[1].grid(which='both')
axes[1].set_title('Vertical Seismic Displ.')

# Reference Kolmogorov scaling
axes[0].loglog(1/freqs, 1.5e0/(freqs_les**(7/3)), 'b--', linewidth=2)

plt.tight_layout()
fig.savefig(os.path.join(les_dir, 'spectra_ref.png'), dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.close(fig)