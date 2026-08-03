#!/usr/bin/env python3

"""
Functions for power spectral density (PSD) calculation
Methods: FFT, Welch, Short-Time FT (STFT), Continuous Wavelet Transform (CWT), fast fCWT package

Author: Qing Ji
"""

# Load python packages

from __future__ import division
import numpy as np
from scipy.signal import convolve2d

import scipy.fft as fft
from scipy.signal import ShortTimeFFT, welch, csd, coherence
from scipy.signal import detrend, savgol_filter
from scipy.interpolate import interp1d
import pycwt as wavelet

import warnings
warnings.filterwarnings('ignore')

# fCWT requires building from source (see fcwt_install/ and top-level README)
# Only the wavelet-based functions (init_fcwt, my_fcwt, my_fxwt, ...) need it
try:
    import fcwt
except ImportError:
    import sys
    fcwt = None
    print("Warning: fCWT not installed. Wavelet functions (init_fcwt, my_fcwt, my_fxwt, ...) "
          "are unavailable. See the fCWT section of the top-level README.", file=sys.stderr)


### Fourier analysis ###

# Hanning window
def hann_taper(data):
    # Hanning window
    return data * np.hanning(len(data)).astype(data.dtype)

# FFT
def my_fft(data, dt, mode='psd', nfft=None):

    # Time series paramters
    if nfft is None:
        Nt = len(data)
        nfft = fft.next_fast_len(Nt)

    # Preprocessing
    _data = detrend(data, type='linear')
    _data = hann_taper(_data)

    # Compensation factor for Hann window
    window_norm = (hann_taper(np.ones(nfft, _data.dtype))**2).sum()

    # Power spectral density
    if mode == 'psd':

        # Positive frequency samples
        freqs = fft.rfftfreq(nfft, d=dt)[1:]

        # Two-sided PSD (not doubled, with energy compensated for window function)
        psd = np.abs(fft.rfft(_data, n=nfft)[1:])**2 * dt / window_norm

        return psd, freqs

    # Fourier spectrum
    elif mode == 'spec':

        # Positive frequency samples
        freqs = fft.rfftfreq(nfft, d=dt)[1:]

        # Spectral amplitude
        spec = fft.rfft(_data, n=nfft)[1:] * dt

        return spec, freqs


# Welch method
def my_welch(data, dt, nperseg, noverlap=None, ave_method='mean', detrend='linear'):

    # Power spectral density
    freqs, psd = welch(data, 1/dt, window='hann', average=ave_method, 
                       nperseg=nperseg, noverlap=noverlap, nfft=None,
                       detrend=detrend, return_onesided=True, scaling='density')

    # Two-sided PSD for postive frequency samples
    freqs = freqs[1:]
    psd = psd[1:] / 2

    return psd, freqs


def my_welch_csd(u1, u2, dt, nperseg, noverlap=None, ave_method='mean'):

    # Welch analysis
    freqs, c12 = csd(u1, u2, 1/dt, window='hann', average=ave_method,
                     nperseg=nperseg, noverlap=noverlap, nfft=None,
                     detrend='linear', return_onesided=True, scaling='density')

    # Two-sided CSD for postive frequency samples
    freqs = freqs[1:]
    c12 = c12[1:] / 2

    return c12, freqs


def my_welch_coh(u1, u2, dt, nperseg, noverlap=None):

    # Welch analysis
    freqs, coh = coherence(u1, u2, 1/dt, window='hann', nperseg=nperseg,
                           noverlap=noverlap, nfft=None, detrend='linear')
    freqs = freqs[1:]
    coh = coh[1:]

    return coh, freqs


# STFT
def my_stft(data, dt, nperseg, hop=1, mode='psd', mfft=None, window=('hann',)):

    # Time series paramters
    Nt = len(data)

    # Power spectral density
    if mode == 'psd':
        # Create STFT class
        SFT = ShortTimeFFT.from_window(window, 1/dt, nperseg, noverlap=nperseg-hop,
                                       scale_to='psd', fft_mode='onesided', mfft=mfft)

        # STFT spectrogram: Two-sided PSD
        psd = SFT.spectrogram(data)[1:, :]
        freqs = SFT.f[1:]
        t_psd = SFT.t(Nt)

        # Remove edge (padding) effect
        ind1 = SFT.lower_border_end[1] - SFT.p_min
        ind2 = SFT.upper_border_begin(Nt)[1] - SFT.p_min
        psd = psd[:, ind1:ind2]
        t_psd = t_psd[ind1:ind2]

        return psd, freqs, t_psd

    # Fourier spectrum
    elif mode == 'spec':
        # Create STFT class
        SFT = ShortTimeFFT.from_window(window, 1/dt, nperseg, noverlap=nperseg-hop,
                                       scale_to=None, fft_mode='onesided', mfft=mfft)

        # STFT spectrum: Not-doubled spectral amplitude
        spec = SFT.stft(data) * dt
        freqs = SFT.f
        t_spec = SFT.t(Nt)

        # Remove edge (padding) effect
        ind1 = SFT.lower_border_end[1] - SFT.p_min
        ind2 = SFT.upper_border_begin(Nt)[1] - SFT.p_min
        spec = spec[:, ind1:ind2]
        t_spec = t_spec[ind1:ind2]

        return spec, freqs, t_spec



### Wavelet analysis: PyCWT ###

# Empirical factor for COI determination
coi_coeff = 2 * 1.6 / np.sqrt(2)

# 2D convolution of matrices x and y
def conv2(x, y, mode='same'):
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)

# Exponents p for the smallest powers of two that satisfy the relation: 2**p >= abs(x)
def nextpow2(x):
    res = np.ceil(np.log2(x))
    return res.astype('int')

# Smoothing function
def smoothCFS(cfs, scales, dt, ns, nt):
    N = cfs.shape[1]
    npad = 2 ** nextpow2(N)
    omega = np.arange(1, np.fix(npad / 2) + 1, 1).tolist()
    omega = np.array(omega) * ((2 * np.pi) / npad)
    omega_save = -omega[int(np.fix((npad - 1) / 2)) - 1:0:-1]
    omega_2 = np.concatenate((0., omega), axis=None)
    omega_2 = np.concatenate((omega_2, omega_save), axis=None)
    omega = np.concatenate((omega_2, -omega[0]), axis=None)
    # Normalize scales by DT because we are not including DT in the angular frequencies here.
    # The smoothing is done by multiplication in the Fourier domain.
    normscales = scales / dt

    # Gaussian window: Standard dev. in time is T * sqrt(2*nt)
    for kk in range(0, cfs.shape[0]):
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk - 1], npad))
        cfs[kk - 1] = smooth[0:N]
    # Convolve the coefficients with a moving average smoothing filter across scales.
    H = 1 / ns * np.ones((ns, 1))

    cfs = conv2(cfs, H)
    return cfs

# Wavelet spectrum / PSD
def my_cwt(u, dt, mode='psd', vpo=12, freqmin=0.02, freqmax=0.1, nptsfreq=200,
           sampling='linear', smooth=False, ns=3, nt=0.25, w0=6.0):
    # Choosing a Morlet wavelet with a central frequency w0 = 6
    mother = wavelet.Morlet(w0)
    # nx represent the number of element in the trace_current array
    nx = np.size(u)
    x = np.transpose(u)
    # Spacing between discrete scales, the default value is 1/12
    dj = 1 / vpo
    # Number of scales less one, -1 refers to the default value which is J = (log2(N * dt / s0)) / dj.
    J = -1
    # Smallest scale of the wavelet, default value is 2*dt
    s0 = 2 * dt  # Smallest scale of the wavelet, default value is 2*dt

    # Creation of the frequency vector that we will use in the continuous wavelet transform
    if sampling == 'linear':
        freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)
    elif sampling == 'log':
        freqlim = np.linspace(np.log10(freqmax), np.log10(freqmin), num=nptsfreq, endpoint=True)
        freqlim = np.power(10, freqlim)
    else:
        freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)

    # Preprocessing data
    x = detrend(x, type='linear')

    # Calculation wavelet transform
    # scales are calculated using the wavelet Fourier wavelength
    # fft : Normalized fast Fourier transform of the input trace
    # fftfreqs : Fourier frequencies for the calculated FFT spectrum.
    cwt, scales, freqs, coi, _, _ = wavelet.cwt(x, dt, dj, s0, J, mother, freqs=freqlim)

    # Modify factors to agree with benchmark on longer trace
    coi = coi / coi_coeff

    # Power spectral density
    if mode == 'psd':
        if smooth:
            scales = np.array([[kk] for kk in scales])
            invscales = np.kron(np.ones((1, nx)), 1 / scales)
            cfs = smoothCFS(invscales * abs(cwt) ** 2, scales, dt, ns, nt)
            cfs = cfs * dt * scales
            coi = coi / np.sqrt(2*nt)   # COI due to temporal smoothing
        else:
            cfs = abs(cwt) ** 2 * dt
        return cfs, freqs, coi

    # Input amplitude
    elif mode == 'spec':
        cwt = cwt / np.sqrt(scales[:, np.newaxis] * dt)
        return cwt, freqs, coi

"""
    elif mode == 'spec':
        if smooth:
            scales = np.array([[kk] for kk in scales])
            invscales = np.kron(np.ones((1, nx)), 1 / scales)
            cwt = smoothCFS(invscales * cwt, scales, dt, ns, nt)
            cwt = cwt * scales
        else:
            if rectify:
                cwt = cwt / np.sqrt(scales[:, np.newaxis]/dt)
            else:
                cwt = cwt * dt
        return cwt, freqs, coi
"""


# Wavelet CSD and Coherence
def my_xwt(u1, u2, dt, vpo=12, freqmin=0.02, freqmax=0.1, nptsfreq=200,
           sampling='linear', smooth=False, ns=3, nt=0.25, w0=6.0):
    # Choosing a Morlet wavelet with a central frequency w0 = 6
    mother = wavelet.Morlet(w0)
    # nx represent the number of element in the trace_current array
    nx = np.min([np.size(u1), np.size(u2)])
    x1 = np.transpose(u1)
    x2 = np.transpose(u2)
    # Spacing between discrete scales, the default value is 1/12
    dj = 1 / vpo
    # Number of scales less one, -1 refers to the default value which is J = (log2(N * dt / so)) / dj.
    J = -1
    # Smallest scale of the wavelet, default value is 2*dt
    s0 = 2 * dt  # Smallest scale of the wavelet, default value is 2*dt

    # Creation of the frequency vector that we will use in the continuous wavelet transform
    if sampling == 'linear':
        freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)
    elif sampling == 'log':
        freqlim = np.linspace(np.log10(freqmax), np.log10(freqmin), num=nptsfreq, endpoint=True)
        freqlim = np.power(10, freqlim)
    else:
        freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)

    # Preprocessing data
    x1 = detrend(x1, type='linear')
    x2 = detrend(x2, type='linear')

    # Calculation of the two wavelet transform independently
    # scales are calculated using the wavelet Fourier wavelength
    # fft : Normalized fast Fourier transform of the input trace
    # fftfreqs : Fourier frequencies for the calculated FFT spectrum.
    cwt1, scales, freqs, coi, _, _ = wavelet.cwt(x1, dt, dj, s0, J, mother, freqs=freqlim)
    cwt2, _, _, _, _, _ = wavelet.cwt(x2, dt, dj, s0, J, mother, freqs=freqlim)

    # Modify factors to agree with benchmark on longer trace
    coi = coi / coi_coeff

    # Cross-wavelet transform operation
    crossCFS = np.conj(cwt1) * cwt2

    # Smoothing parameters
    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)

    if smooth:
        cfs1 = smoothCFS(invscales * abs(cwt1) ** 2, scales, dt, ns, nt)
        cfs2 = smoothCFS(invscales * abs(cwt2) ** 2, scales, dt, ns, nt)
        crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
        Wcoh = abs(crossCFS) ** 2 / (cfs1 * cfs2)
        coi = coi / np.sqrt(2*nt)   # COI due to temporal smoothing

        # Consistent with Fourier analysis
        cfs1 = cfs1 * dt * scales
        cfs2 = cfs2 * dt * scales
        crossCFS = crossCFS * dt * scales

    else:
        # Consistent with Fourier analysis
        cfs1 = abs(cwt1) ** 2 * dt
        cfs2 = abs(cwt2) ** 2 * dt
        crossCFS = crossCFS * dt

        crossCFS_sm = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
        cfs1_sm = smoothCFS(invscales * cfs1, scales, dt, ns, nt)
        cfs2_sm = smoothCFS(invscales * cfs2, scales, dt, ns, nt)
        Wcoh = abs(crossCFS_sm) ** 2 / (cfs1_sm * cfs2_sm)

    return cfs1, cfs2, crossCFS, Wcoh, freqs, coi


### Fast fCWT (Arts & Van Den Broek, 2022, Nat Comput Sci) ###

# Initialize fCWT object
def init_fcwt(w0=6.0, Nt=16384, nthreads=4, fftw_plan='FFTW_MEASURE'):

    # Mother wavelet (Morlet only)
    mother = fcwt.Morlet(w0/(2*np.pi))

    # Initialize fCWT object
    use_optimization_plan = True
    use_normalization = True
    fcwt_obj = fcwt.FCWT(mother, nthreads, use_optimization_plan, use_normalization)
    fcwt_obj.create_FFT_optimization_plan(Nt, fftw_plan)

    return fcwt_obj, mother


# Wavelet spectrum / PSD
def my_fcwt(fcwt_obj, mother, u, dt, mode='psd', freqmin=0.02, freqmax=0.1, nptsfreq=200,
            sampling='linear', smooth=False, ns=3, nt=0.25):

    # Prepare for fCWT
    Nt = len(u)                 # Length of signal
    fs = 1 / dt                 # Sampling rate (Hz)
    w0 = mother.fb * (2*np.pi)  # Morlet wavelet parameter
    flambda = (4*np.pi) / (w0 + np.sqrt(2 + w0 ** 2))   # Fourier wavelength
    f_factor = 0.978            # Factor to account for frequency (currently)
    data = u.astype('single')

    # Evenly-spaced samples in linear or log scales
    if sampling == 'linear':
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LINFREQS, fs, freqmin, freqmax, nptsfreq)
    elif sampling == 'log':
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LOGSCALES, fs, freqmin, freqmax, nptsfreq)
    else:
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LINFREQS, fs, freqmin, freqmax, nptsfreq)

    # Frequency and scale samples
    freqs = np.zeros((nptsfreq), dtype='single')
    scale_obj.getFrequencies(freqs)
    freqs = freqs / f_factor
    scales = 1 / (flambda * freqs)

    # Preprocessing data
    data = detrend(data, type='linear')

    # Calculate continuous wavelet transform
    cwt = np.zeros((nptsfreq, data.size), dtype=np.complex64)
    fcwt_obj.cwt(data, scale_obj, cwt)

    # Determines the cone-of-influence (COI). Note it is returned as a function
    # of time in Fourier periods. Uses triangular Bartlett window with
    # non-zero end-points
    # Modify factors to agree with benchmark on longer trace
    coi = (Nt / 2 - np.abs(np.arange(0, Nt) - (Nt - 1) / 2))
    coi = flambda / np.sqrt(2) * dt * coi / coi_coeff

    # Power spectral density
    if mode == 'psd':
        if smooth:
            _scales = np.array([[kk] for kk in scales])
            cfs = smoothCFS(abs(cwt)**2, _scales, dt, ns, nt)
            cfs = cfs * scales[:, np.newaxis]
            coi = coi / np.sqrt(2*nt)   # COI due to temporal smoothing
        else:
            cfs = abs(cwt)**2 * scales[:, np.newaxis]
        return cfs, freqs, coi

    # Input amplitude
    elif mode == 'spec':
        cwt = cwt / dt
        return cwt, freqs, coi


# Wavelet CSD and Coherence
def my_fxwt(fcwt_obj, mother, u1, u2, dt, freqmin=0.02, freqmax=0.1, nptsfreq=200,
            sampling='linear', smooth=False, ns=3, nt=0.25):

    # Prepare for fCWT
    Nt = len(u1)                # Length of signal
    if len(u2) != Nt:
        raise ValueError('Input signals not having the same length!')
    fs = 1 / dt                 # Sampling rate (Hz)
    w0 = mother.fb * (2*np.pi)  # Morlet wavelet parameter
    flambda = (4*np.pi) / (w0 + np.sqrt(2 + w0 ** 2))   # Fourier wavelength
    f_factor = 0.978            # Factor to account for frequency (currently)
    data1 = u1.astype('single')
    data2 = u2.astype('single')

    # Evenly-spaced samples in linear or log scales
    if sampling == 'linear':
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LINFREQS, fs, freqmin, freqmax, nptsfreq)
    elif sampling == 'log':
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LOGSCALES, fs, freqmin, freqmax, nptsfreq)
    else:
        scale_obj = fcwt.Scales(mother, fcwt.FCWT_LINFREQS, fs, freqmin, freqmax, nptsfreq)

    # Frequency and scale samples
    freqs = np.zeros((nptsfreq), dtype='single')
    scale_obj.getFrequencies(freqs)
    freqs = freqs / f_factor
    scales = 1 / (flambda * freqs)

    # Preprocessing data
    data1 = detrend(data1, type='linear')
    data2 = detrend(data2, type='linear')

    # Calculate continuous wavelet transform
    cwt1 = np.zeros((nptsfreq, data1.size), dtype=np.complex64)
    fcwt_obj.cwt(data1, scale_obj, cwt1)
    cwt2 = np.zeros((nptsfreq, data2.size), dtype=np.complex64)
    fcwt_obj.cwt(data2, scale_obj, cwt2)

    # Wavelet PSD & CSD
    cfs1 = abs(cwt1) ** 2
    cfs2 = abs(cwt2) ** 2
    crossCFS = np.conj(cwt1) * cwt2

    # Determines the cone-of-influence (COI). Note it is returned as a function
    # of time in Fourier periods. Uses triangular Bartlett window with
    # non-zero end-points
    # Modify factors to agree with benchmark on longer trace
    coi = (Nt / 2 - np.abs(np.arange(0, Nt) - (Nt - 1) / 2))
    coi = flambda / np.sqrt(2) * dt * coi / coi_coeff

    # Coherence analysis
    _scales = np.array([[kk] for kk in scales])
    _cfs1 = smoothCFS(cfs1, _scales, dt, ns, nt)
    _cfs2 = smoothCFS(cfs2, _scales, dt, ns, nt)
    _crossCFS = smoothCFS(crossCFS, _scales, dt, ns, nt)
    Wcoh = abs(_crossCFS) ** 2 / (_cfs1 * _cfs2)

    if smooth:
        _cfs1 = _cfs1 * scales[:, np.newaxis]
        _cfs2 = _cfs2 * scales[:, np.newaxis]
        _crossCFS = _crossCFS * scales[:, np.newaxis]
        coi = coi / np.sqrt(2*nt)   # COI due to temporal smoothing
        return _cfs1, _cfs2, _crossCFS, Wcoh, freqs, coi

    else:
        cfs1 = cfs1 * scales[:, np.newaxis]
        cfs2 = cfs2 * scales[:, np.newaxis]
        crossCFS = crossCFS * scales[:, np.newaxis]
        return cfs1, cfs2, crossCFS, Wcoh, freqs, coi


def my_fcwt_coh(mother, csd12, cfs1, cfs2, dt, freqs, smooth=True, ns=3, nt=0.25):

    # Obtain scales
    w0 = mother.fb * (2*np.pi)  # Morlet wavelet parameter
    flambda = (4*np.pi) / (w0 + np.sqrt(2 + w0 ** 2))   # Fourier wavelength
    scales = 1 / (flambda * freqs)

    # Coherence analysis
    if smooth:
        _scales = np.array([[kk] for kk in scales])
        _cfs1 = smoothCFS(cfs1, _scales, dt, ns, nt)
        _cfs2 = smoothCFS(cfs2, _scales, dt, ns, nt)
        _crossCFS = smoothCFS(csd12, _scales, dt, ns, nt)
        Wcoh = abs(_crossCFS) ** 2 / (_cfs1 * _cfs2)
    else:
        Wcoh = abs(csd12) ** 2 / (cfs1 * cfs2)

    return Wcoh


### Other functions & Applications ###

# Function: Spectral smoothing
def smooth_spectrum(f_in, spec_in, npts=10, smooth_mode='nearest', polyorder=1,
                    interp=True, f_out=None):

    # Interp. in freq. domain
    if interp:
        func = interp1d(f_in, spec_in, fill_value="extrapolate")
        if f_out is None:
            f_out = 10 ** np.linspace(np.log10(np.min(f_in)), np.log10(np.max(f_in)), len(f_in))
        spec_out = func(f_out)

    # Smoothing
    spec_out = savgol_filter(spec_out, npts, polyorder=polyorder, mode=smooth_mode)
    return spec_out, f_out


# Function: Create period bins
def setup_period_bins(period_range, octave_step=0.125, bin_step=1):
    """
    Setup period bins for smoothing the spectrum
    For meaning of parameters, see McNamara & Buland (2004, BSSA)

    Parameters:
        period_range (list): Period range [T1, T2] with T1 < T2
        octave_step (float): Period step width in log2 scale (default: 0.125)
        bin_step (float): Bin width in log2 scale (default: 1)

    Returns:
        period_samps (3*Np array): Left edge, center and right edge of period bins
    """

    # Period range (check T1 < T2)
    period_range = np.sort(period_range)

    # Create period bins
    period_center = 2 ** np.arange(np.log2(period_range[0]), np.log2(period_range[1])+octave_step, octave_step)
    period_left = period_center / 2**(bin_step/2)
    period_right = period_center * 2**(bin_step/2)

    return np.vstack([period_left, period_center, period_right])


# Function: Find horizontal direction with maximum coherence
def find_max_coh_angle(uz, un, ue, dt, freq_range=[0.01, 0.05], nperseg=1024, dtheta=1):

    # Initialize arrays
    theta_arr = np.arange(0, 180, dtheta)
    coh_arr = np.zeros(theta_arr.shape)

    # Frequency samples
    freqs, coh = coherence(uz, un, fs=1/dt, window='hann', nperseg=int(nperseg/dt),
                           noverlap=int(nperseg/dt*0.75), detrend='linear')
    f_mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
    coh_arr[0] = np.median(coh[f_mask])

    for i in range(1, len(theta_arr)):

        # Rotated horizontal component
        theta = np.deg2rad(theta_arr[i])
        ur = un * np.cos(theta) + ue * np.sin(theta)

        # Coherence
        _, coh = coherence(uz, ur, fs=1/dt, window='hann', nperseg=int(nperseg/dt),
                           noverlap=int(nperseg/dt*0.75), detrend='linear')
        coh_arr[i] = np.median(coh[f_mask])

    return theta_arr, coh_arr
