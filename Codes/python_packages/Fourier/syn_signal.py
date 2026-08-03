#!/usr/bin/env python3

"""
Functions for generating signals to test Fourier analysis methods
Synthetic signals: 
    1. Sinusoidal 
    2. Noise from band spectrum 
    3. Wave packet from band spectrum
    4. Seismic ambient noise 
    5. Synthetic signal for fCWT 
    6. Human EEG dataset
    7. Synthetic signal with dispersion

Author: Qing Ji
"""

# Load python packages

from __future__ import division
import os
import numpy as np
import pandas as pd

import scipy.fft as fft
from scipy.integrate import cumulative_trapezoid

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

from seismic.seis_download import download_trace

import warnings
warnings.filterwarnings('ignore')

# Data directory
data_dir = '/Users/qingji/Documents/Datasets/cwt_syn'


### Generate synthetic signals ###

# Sinusoidal signal
def data_sine(dt, Nt, f_arr, amp_arr, phi_arr=None):
    
    # Time axis [s]
    t = np.arange(0, Nt) * dt

    # Sinusoidal signal
    f_arr = np.array(f_arr).reshape(1, -1)
    amp_arr = np.array(amp_arr).reshape(1, -1)
    if phi_arr is None:
        data = np.sum(amp_arr*np.sin(2*np.pi*f_arr*t[:, None]), axis=1)
    else:
        # Phase delay in time
        data = np.sum(amp_arr*np.sin(2*np.pi*f_arr*t[:, None] - phi_arr), axis=1)
    
    return data, t


# Noise signal from a band spectrum
def data_noise(dt, Nt, fc_arr, amp=1):
    
    # Time axis [s]
    t = np.arange(0, Nt) * dt

    # Frequency samples [Hz, rad/s]
    f = fft.rfftfreq(Nt, dt)
    w = 2*np.pi * f
    df = f[1] - f[0]
    
    # Band spectrum
    f1, f2, f3, f4 = fc_arr[0], fc_arr[1], fc_arr[2], fc_arr[3]
    psd = amp * band_spec(f, f1, f2, f3, f4)

    # Random phase
    phi = np.random.uniform(0, 2*np.pi, len(f))

    # Synthetic noise
    syn_noise = np.sum(2*np.sqrt(psd[1:, None]*df) * np.cos(w[1:, None] @ t[None] + phi[1:, None]), axis=0)
    
    return syn_noise, t, f, psd


# Wave packet from a band spectrum
def data_band(dt, Nt, fc_arr, amp=1):
    
    # Time axis [s]
    t = np.arange(0, Nt) * dt

    # Frequency samples [Hz, rad/s]
    f = fft.rfftfreq(Nt, dt)
    w = 2*np.pi * f
    df = f[1] - f[0]
    
    # Band spectrum
    f1, f2, f3, f4 = fc_arr[0], fc_arr[1], fc_arr[2], fc_arr[3]
    spec = amp * band_spec(f, f1, f2, f3, f4)

    # Synthetic signal
    syn_data = np.sum((2*spec[1:, None]*df) * np.cos(w[1:, None] @ (t[None] - Nt*dt/2)), axis=0)
    
    return syn_data, t, f, spec


# Seismic ambient noise
def data_seis(sta_char=['TA','645A'], channel='LHZ', resp_to='DISP', extra_portion=1/12,
              ts=UTCDateTime(2012,8,29,7), te=UTCDateTime(2012,8,29,8)):
    
    # Seismic data
    trace, _ = download_trace(sta_char, [ts, te], channel, pre_filt='default',
                              extra_portion=extra_portion, client=Client('IRIS'), resp_to=resp_to)
    dt, Nt = trace.stats.delta, trace.stats.npts
    t = np.arange(0, Nt*dt, dt)
    
    return trace, t


# Band spectrum
def band_spec(f, f1, f2, f3, f4, water_level=1e-4):
    
    A = np.zeros(f.shape)
    
    mask1 = (f<f2)
    mask2 = (f>f3)
    mask3 = (f>=f2) & (f<=f3)
    
    # Cosine taper
    # A[mask1] = np.cos(np.pi/(f2-f1)/2 * (f[mask1]-f2)) ** 2
    # A[mask2] = np.cos(np.pi/(f4-f3)/2 * (f[mask2]-f3)) ** 2
    # A[mask3] = 1
    
    # Gaussian taper
    A[mask1] = np.exp(-0.5 * ((f[mask1] - f2) / ((f2 - f1) / 4))**2)
    A[mask2] = np.exp(-0.5 * ((f[mask2] - f3) / ((f4 - f3) / 4))**2)
    A[mask3] = 1
    
    # Water level
    min_level = water_level
    A = (A * (1 - min_level)) + min_level
    
    return A


# fCWT synthetic signal
def data_fcwt(noise=False):

    # Read csv file
    csvFile = os.path.join(data_dir, 'synthetic_dataset.csv')
    df = pd.read_csv(csvFile)

    # Read synthetic data
    t = df.iloc[:, 0].values
    if noise:
        data = df.iloc[:, 2].values
    else:
        data = df.iloc[:, 1].values

    return data, t


# Human EEG dataset
def data_EEG():

    # Read csv file
    csvFile = os.path.join(data_dir, 'EEG_dataset.csv')
    df = pd.read_csv(csvFile)

    # Read EEG dataset
    t = df.iloc[:, 0].values
    data = df.iloc[:, 1].values

    return data, t


### Dispersive wave ###

# Phase -> Group
def vph_to_vg(omega, vph):
    """Covert vph into vg using finite differences."""
    if callable(vph):
        vph_ = np.asarray(vph(omega))
    else:
        vph_ = np.asarray(vph)
    k = omega / vph_
    return np.gradient(omega, k)

# Group -> Phase
def vg_to_vph(omega, vg, vp0):
    """Covert vg into vph using finite differences."""
    if callable(vg):
        vg_ = np.asarray(vg(omega))
    else:
        vg_ = np.asarray(vg)
    k = np.concatenate(([omega[0]/vp0], cumulative_trapezoid(1.0/vg_, omega)))
    return omega / k

# Synthetic dispersive wave
def dispersive_wave(t, x, omega, vph, amp=None):
    """Assemble dispersive wave; *vph* can be array or callable.
    The signal is returned on the time grid *t*.
    """
    if callable(vph):
        vph_ = np.asarray(vph(omega))
    else:
        vph_ = np.asarray(vph)

    if amp is None:
        amp = np.ones_like(omega)
    if amp.shape != omega.shape:
        raise ValueError("Input amplitude must match omega shape.")
    
    k = omega / vph_
    exp_term = np.exp(1j * (k[:, None] * x - omega[:, None] * t[None, :]))
    wave = (amp[:, None] * exp_term).sum(axis=0)
    return np.real(wave)
