#!/usr/bin/env python3

"""
Long period pressure PSD for misfit

Author: Qing Ji
"""

# Load python packages
import os
import numpy as np
import pandas as pd


# Directory for CM1 modeling outputs
les_dir = '/path/to/simulations'    # MODIFY THIS PATH

# Get all simulation directories
sim_dirs = sorted([d for d in os.listdir(les_dir) 
			    if os.path.isdir(os.path.join(les_dir, d)) and d.endswith("dense")])

# Store pressure PSD
p_psd_arr = []

for sim_dir in sim_dirs:
    print(f"Processing directory: {sim_dir}")
    
    # Read spectra
    psd_les = np.load(os.path.join(sim_dir, 'spectra_height', f'spectra_prs_15m.npz'))
    
    # Spectral amplitude
    mask = (1/psd_les['freqs'] >= 30) & (1/psd_les['freqs'] <= 200)
    p_psd = np.mean(psd_les['psd'][0][mask])
    print(f'Pressure PSD: {p_psd:.2g} Pa^2/Hz')
    p_psd_arr.append(p_psd)
    
# Save to CSV file
df = pd.DataFrame({"case": sim_dirs, "p_psd": p_psd_arr})
df.to_csv("psd_30-200s.csv", index=False)