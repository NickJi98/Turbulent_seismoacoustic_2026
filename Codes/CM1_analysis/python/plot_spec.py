#!/usr/bin/env python3

"""
Plot CM1 spectra at different heights compared with observation

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


# Observation
obs_spec = np.load('/path/to/Data/Isaac/obs_spectra.npz')
ind_t = np.where(obs_spec['h_tower'] == 15)[0][0]
freqs = obs_spec['f_tower']

# Directory for CM1 spectra
les_dir = './'
spec_dir = 'spectra_height'

### Wind spectra ###
height_list = [5, 15, 25, 35, 45, 95, 195]
print("Heights:", height_list)

# Plot wind spectra
fig, axes = plt.subplots(1, 2, figsize=(5, 2.5), sharex=True)
for height in height_list:
    try:
        psd_les = np.load(os.path.join(les_dir, spec_dir, f'spectra_wind_{height}m.npz'))
    except:
        continue
    axes[0].loglog(1/psd_les['freqs'], psd_les['psd_u1'][0], '-', label=f'{height} m', lw=1)
    axes[1].loglog(1/psd_les['freqs'], psd_les['psd_u2'][0], '-', label=f'{height} m', lw=1)

# Wind tower
axes[0].loglog(1/freqs, obs_spec['Suu'][ind_t], 'k-', label='Tower 15 m', lw=2)
axes[1].loglog(1/freqs, obs_spec['Svv'][ind_t], 'k-', label='Tower 15 m', lw=2)

# Modification
axes[0].set_xlabel('Period (s)')
axes[0].set_ylabel('PSD ((m/s)$^2$/Hz)')
axes[0].set_xlim(1/np.flip([5e-3, 0.25]))
axes[0].set_ylim([1e0, 1e3])
axes[0].grid(which='both')
axes[0].legend(loc='lower right', prop={'size': 4})
axes[0].set_title('Along-wind', fontsize=10)
axes[1].set_xlabel('Period (s)')
axes[1].set_ylabel('PSD ((m/s)$^2$/Hz)')
axes[1].set_ylim([5e-1, 5e2])
axes[1].grid(which='both')
axes[1].set_title('Cross-wind', fontsize=10)

# Reference Kolmogorov scaling
axes[0].loglog(1/freqs, 2e-1/(freqs**(5/3)), 'b--', linewidth=2)
axes[1].loglog(1/freqs, 2e-1/(freqs**(5/3)), 'b--', linewidth=2)

plt.tight_layout()
fig.savefig(os.path.join(les_dir, spec_dir, 'spectra_wind.png'), dpi=200, 
            bbox_inches='tight', pad_inches=0.1)
plt.close(fig)


### Temperature & pressure spectra ###
fig, axes = plt.subplots(1, 2, figsize=(5, 2.5), sharex=True)
for height in height_list:
    psd_les = np.load(os.path.join(les_dir, spec_dir, f'spectra_th_{height}m.npz'))
    axes[0].loglog(1/psd_les['freqs'], psd_les['psd'][0], '-', label=f'{height} m', lw=1)
    psd_les = np.load(os.path.join(les_dir, spec_dir, f'spectra_prs_{height}m.npz'))
    axes[1].loglog(1/psd_les['freqs'], psd_les['psd'][0], '-', label=f'{height} m', lw=1)

# Observation
axes[0].loglog(1/freqs, obs_spec['Stt'][ind_t], 'k-', label='Tower 15 m', lw=2)
axes[1].loglog(1/obs_spec['f_seis'], obs_spec['Spp'], 'k-', label='Infrasound', lw=2)

# Modification
axes[0].set_xlabel('Period (s)')
axes[0].set_ylabel('PSD (K$^2$/Hz)')
axes[0].set_xlim(1/np.flip([5e-3, 0.25]))
axes[0].set_ylim([2e-3, 2e0])
axes[0].grid(which='both')
axes[0].legend(loc='lower right', prop={'size': 4})
axes[0].set_title('Temperature', fontsize=10)
axes[1].set_xlabel('Period (s)')
axes[1].set_ylabel('PSD (Pa$^2$/Hz)')
axes[1].set_ylim([1e1, 2e3])
axes[1].grid(which='both')
axes[1].set_title('Pressure', fontsize=10)

# Reference Kolmogorov scaling
axes[0].loglog(1/freqs, 1e-3/(freqs**(5/3)), 'b--', linewidth=2)
axes[1].loglog(1/freqs, 8e-1/(freqs**(7/3)), 'b--', linewidth=2)

plt.tight_layout()
fig.savefig(os.path.join(les_dir, spec_dir, 'spectra_tp.png'), dpi=200, 
            bbox_inches='tight', pad_inches=0.1)
plt.close(fig)