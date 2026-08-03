#!/usr/bin/env python3

"""
Gather CM1 spectra at different heights to a single file

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np


# Directory for CM1 spectra
les_dir = './V42_Cd22_dense'
spec_dir = 'spectra_height'

# List of height levels
height_list = [5, 15, 25, 35, 45, 95, 195]
print("Heights:", height_list)

# Gather wind spectra
psd_u1_list, psd_u2_list = [], []
freqs = None

for h in height_list:
    fn = os.path.join(les_dir, spec_dir, f'spectra_wind_{h}m.npz')
    if not os.path.exists(fn):
        print(f"Skipping {fn} (missing)")
        continue
    
    data = np.load(fn)
    if freqs is None:
        freqs = data['freqs']
    
    psd_u1_list.append(data['psd_u1'])
    psd_u2_list.append(data['psd_u2'])

# Stack into (nz, 3, nfreqs)
psd_u1 = np.stack(psd_u1_list, axis=0)
psd_u2 = np.stack(psd_u2_list, axis=0)

# Save wind spectra
out_fn = os.path.join(les_dir, spec_dir, "spectra_wind_all.npz")
np.savez(out_fn, h=np.array(height_list), freqs=freqs,
         psd_u1=psd_u1, psd_u2=psd_u2)


# Gather temperature spectra
psd_th_list, freqs = [], None

for h in height_list:
    fn = os.path.join(les_dir, spec_dir, f'spectra_th_{h}m.npz')
    if not os.path.exists(fn):
        print(f"Skipping {fn} (missing)")
        continue
    
    data = np.load(fn)
    if freqs is None:
        freqs = data['freqs']
    psd_th_list.append(data['psd'])
    
# Save temperature spectra
psd_th = np.stack(psd_th_list, axis=0)
out_fn = os.path.join(les_dir, spec_dir, "spectra_th_all.npz")
np.savez(out_fn, h=np.array(height_list), freqs=freqs, psd_th=psd_th)


# Gather pressure spectra
psd_prs_list, freqs = [], None

for h in height_list:
    fn = os.path.join(les_dir, spec_dir, f'spectra_prs_{h}m.npz')
    if not os.path.exists(fn):
        print(f"Skipping {fn} (missing)")
        continue
    
    data = np.load(fn)
    if freqs is None:
        freqs = data['freqs']
    psd_prs_list.append(data['psd'])

# Save pressure spectra
psd_prs = np.stack(psd_prs_list, axis=0)
out_fn = os.path.join(les_dir, spec_dir, "spectra_prs_all.npz")
np.savez(out_fn, h=np.array(height_list), freqs=freqs, psd_prs=psd_prs)
