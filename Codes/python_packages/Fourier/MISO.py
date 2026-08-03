#!/usr/bin/env python3

"""
Multi-Input Single-Output (MI/SO) analysis

Author: Qing Ji
"""

# Load python packages

from __future__ import division
import numpy as np
from scipy.signal import convolve2d

import scipy.fft as fft
from scipy.signal import ShortTimeFFT, welch, csd, coherence

from . import spectrum as myspec


### MISO analysis ###

# Two component
def miso_2input(psd, csd):

    # Input: 1, 2 (decreasing coherence)
    # Output: o
    
    # Linear estimators
    L1o = csd['1o'] / psd['1']
    L12 = csd['12'] / psd['1']
    C2o_1 = csd['2o'] - np.conj(L12) * csd['1o']
    C22_1 = psd['2'] - abs(csd['12'])**2 / psd['1']
    L2o = C2o_1 / C22_1
    C1o_2 = csd['1o'] - csd['12'] / psd['2'] * csd['2o']
    C11_2 = psd['1'] - abs(csd['12'])**2 / psd['2']

    # Transfer function
    H1o = L1o - L12 * L2o

    # Partial coherence
    coh_2o_1 = np.abs(C2o_1)**2 / (C22_1 * psd['o'])
    coh_1o_2 = np.abs(C1o_2)**2 / (C11_2 * psd['o'])

    # Output
    tf = {'1o': H1o, 'L1o': L1o, '12': L12, '2o': L2o}
    par_coh = {'2o_1': coh_2o_1, '1o_2': coh_1o_2}

    return tf, par_coh


# Three component
def miso_3input(psd, csd):

    # Input: 1, 2, 3 (decreasing coherence)
    # Output: o
    
    # Linear estimators
    L1o = csd['1o'] / psd['1']
    L12 = csd['12'] / psd['1']
    L13 = csd['13'] / psd['1']
    C2o_1 = csd['2o'] - np.conj(L12) * csd['1o']
    C22_1 = psd['2'] - abs(csd['12'])**2 / psd['1']
    C23_1 = csd['23'] - np.conj(L12) * csd['13']
    C3o_1 = csd['3o'] - np.conj(L13) * csd['1o']
    C33_1 = psd['3'] - abs(csd['13'])**2 / psd['1']
    L2o = C2o_1 / C22_1
    L23 = C23_1 / C22_1
    C3o_21 = C3o_1 - np.conj(L23) * C2o_1
    C33_21 = C33_1 - abs(C23_1)**2 / C22_1
    L3o = C3o_21 / C33_21

    # Transfer function
    H2o = L2o - L23 * L3o
    H1o = L1o - L12 * H2o - L13 * L3o

    # Partial coherence
    coh_2o_1 = np.abs(C2o_1)**2 / (C22_1 * psd['o'])
    coh_3o_21 = np.abs(C3o_21)**2 / (C33_21 * psd['o'])

    # Output
    tf = {'1o': H1o, '2o': H2o, '3o': L3o, '12': L12, '13': L13, 
          '23': L23, 'L1o': L1o, 'L2o': L2o}
    par_coh = {'2o_1': coh_2o_1, '3o_21': coh_3o_21}

    return tf, par_coh