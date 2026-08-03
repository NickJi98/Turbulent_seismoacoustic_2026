#!/usr/bin/env python3

"""
PDF Class for PSDPDF calculation

Author: Qing Ji
"""

# Load python packages
import fnmatch
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from obspy import UTCDateTime

from scipy.signal import welch
from . import spectrum as myspec

# fCWT requires building from source (see fcwt_install/ and top-level README)
try:
    import fcwt
except ImportError:
    import sys
    fcwt = None
    print("Warning: fCWT not installed. PSDPDF wavelet methods are unavailable. "
          "See the fCWT section of the top-level README.", file=sys.stderr)


### PDF Class ###

class myPDF:

    # Variables to save in .npz output file
    NPZ_STORE_KEYS_BASIC = ['dt', 'window', 'method', 'nfreqs', 'channel', 'resp_comp',
                            'win_shift', 'period_bins', 'var_bins', 
                            'scale_T', 'scale_var', 'histogram']
    NPZ_STORE_KEYS_PSD = ['nperseg', 'w0']
    NPZ_STORE_KEYS_STR_TYPES = ['method', 'scale_T', 'scale_var', 'channel', 'resp_comp']
    NPZ_STORE_KEYS_FIX = ['histogram', 'nperseg', 'w0']

    # A mapping of values for storing info in the NPZ file. This is needed
    # because some types are not loadable without allowing pickle.
    NPZ_SIMPLE_TYPE_MAP = {None: ''}
    NPZ_SIMPLE_TYPE_MAP_R = {v: i for i, v in NPZ_SIMPLE_TYPE_MAP.items()}


    def __init__(self, scale_T='log', scale_var='linear'):
        """
        Initialize the myPDF class.

        Parameters:
        - scale (str): The scale of axes, 'linear' or 'log'.
        """
        self.dt = -1
        self.fs = -1
        self.window = -1
        self.nwin = -1
        self.win_shift = -1
        self.method = 'welch'
        self.nfreqs = -1
        self.channel = 'XXX'
        self.resp_comp = 'ACC'

        self.period_bins = -1
        self.period_bins_centers = -1
        self.var_bins = -1
        self.var_bins_centers = -1
        self.histogram = None
        self.ncount = 0
        self.scale_T = scale_T  # 'linear' or 'log'
        self.scale_var = scale_var  # 'linear' (dB) or 'log' (X^2/Hz)

        # Welch parameters
        self.nperseg = None

        # fCWT parameters
        self.w0 = None

    def set_data_stats(self, channel, window, resp_comp='ACC'):
        if fnmatch.fnmatch(channel, 'L*'):
            self.fs = 1
        elif fnmatch.fnmatch(channel, 'B*'):
            self.fs = 40
        else:
            raise ValueError(f"Band code {channel[0]} not implemented. To be added.")

        self.channel = channel
        self.dt = 1 / self.fs
        self.window = window
        self.nwin = self.window * self.fs

        if fnmatch.fnmatch(channel, '*H[EN12Z]'):
            self.resp_comp = resp_comp
        elif fnmatch.fnmatch(channel, '*D[FH]'):
            self.resp_comp = 'DEF'
        else:
            self.resp_comp = 'DEF'

    def set_period_bins(self, period_samps):
        """
        Set the period bins.

        Parameters:
        - period_samps (array): An array of period bin center values.
        """
        self.period_bins_centers = np.sort(period_samps)
        self.period_bins = self.convert_to_bin_edges(period_samps, self.scale_T)
        self.nfreqs = len(period_samps)

    def set_var_bins(self, var_bins):
        """
        Set the variable bins.

        Parameters:
        - var_bins (array): An array of variable bin edge values.
        """
        self.var_bins = np.array(var_bins)
        self.var_bins_centers = self.convert_to_bin_centers(var_bins, self.scale_var)
        
    def convert_to_bin_edges(self, center_values, scale):
        """
        Convert the center values of bins to bin edges in log10 scale or linear scale.

        Parameters:
        - center_values (array): Array containing the center values of each bin.
        - scale (str): The scale for bin conversion, 'linear' or 'log'.

        Returns:
        - bin_edges (array): Array containing the bin edges in log10 scale or linear scale.
        """
        if scale == 'log':
            bin_step = center_values[1] / center_values[0]
            bin_edges = center_values / np.sqrt(bin_step)
            bin_edges = np.append(bin_edges, bin_edges[-1]*bin_step)
            
        elif scale == 'linear':
            bin_step = center_values[1] - center_values[0]
            bin_edges = center_values - bin_step / 2
            bin_edges = np.append(bin_edges, bin_edges[-1]+bin_step)
            
        else:
            raise ValueError("Invalid scale. Use 'linear' or 'log'.")

        return bin_edges
    
    def convert_to_bin_centers(self, edge_values, scale):
        """
        Convert the bin edges to bin centers in log10 scale or linear scale.

        Parameters:
        - edge_values (array): Array containing the bin edge values.
        - scale (str): The scale for bin conversion, 'linear' or 'log'.

        Returns:
        - bin_centers (array): Array containing the bin centers in log10 scale or linear scale.
        """
        if scale == 'log':
            bin_step = edge_values[1] / edge_values[0]
            bin_centers = edge_values[1:] / np.sqrt(bin_step)
            
        elif scale == 'linear':
            bin_step = edge_values[1] - edge_values[0]
            bin_centers = edge_values[1:] - bin_step / 2
            
        else:
            raise ValueError("Invalid scale. Use 'linear' or 'log'.")

        return bin_centers
    
    @staticmethod
    def sliding_window(data, window_size, step_size):
        """
        Generator for sliding window.

        Parameters:
        - data (array): Input array to slide over.
        - window_size (int): The size of the window.
        - step_size (int): The step size for sliding.

        Yields:
        - array-like: The next window of data.
        """
        for start in range(0, len(data) - window_size + 1, step_size):
            yield data[start:start + window_size]

    def get_response(self, resp_obj):
        self.resp = np.abs(resp_obj.response.get_evalresp_response_for_frequencies(
            1/self.period_bins_centers, output=self.resp_comp)) ** 2 
    
    def init_fft_analyzer(self, method='welch', win_shift=3600, nperseg=None, w0=6.0):

        if method == 'welch':
            self.method = 'welch'
            self.win_shift = win_shift
            self.nperseg = nperseg

            # Frequency samples from PSD calculation
            freqs, _ = welch(np.ones(self.nwin,), fs=self.fs, nperseg=self.nperseg, 
                             noverlap=int(0.75*self.nperseg))
            self.psd_freqs = freqs[1:]

        elif method == 'fcwt':
            self.method = 'fcwt'
            self.win_shift = win_shift
            self.w0 = w0
            fcwt_param = {'w0': self.w0, 'Nt': self.nwin, 
                          'nthreads': 4, 'fftw_plan': 'FFTW_MEASURE'}
            self.fcwt_obj, self.morlet = myspec.init_fcwt(**fcwt_param)
            
            # Frequency samples from PSD calculation
            self.Tmin, self.Tmax = self.period_bins_centers[0], self.period_bins_centers[-1]
            if self.scale_T == 'linear':
                scale_obj = fcwt.Scales(self.morlet, fcwt.FCWT_LINFREQS, self.fs, 
                                        1/self.Tmax, 1/self.Tmin, self.nfreqs)
            elif self.scale_T == 'log':
                scale_obj = fcwt.Scales(self.morlet, fcwt.FCWT_LOGSCALES, self.fs, 
                                        1/self.Tmax, 1/self.Tmin, self.nfreqs)
            freqs = np.zeros((self.nfreqs), dtype='single')
            scale_obj.getFrequencies(freqs)
            self.psd_freqs = freqs

        else:
            raise ValueError("Method should be either 'welch' or 'fcwt'!")
    
    def calc_psd(self, signal):

        # Welch method
        if self.method == 'welch':
            iter_obj = self.sliding_window(signal, self.nwin, self.win_shift*self.fs)
            i, Nseg = 0, len(list(iter_obj))
            psd = np.full((Nseg, self.nfreqs), np.nan)

            iter_obj = self.sliding_window(signal, self.nwin, self.win_shift*self.fs)
            for segment in iter_obj:
                _psd, _ = myspec.my_welch(segment, self.dt, self.nperseg, 
                                          int(0.75*self.nperseg), ave_method='mean')
                
                # Binning. Use geometric mean for averaging in log scale
                for j, per_left, per_right in zip(range(self.nfreqs), self.period_bins[:-1], 
                                                  self.period_bins[1:]):
                    ind = (per_left <= 1/self.psd_freqs) & (1/self.psd_freqs <= per_right)
                    psd[i, j] = stats.gmean(_psd[ind])
                i += 1

            # Remove instrument response
            psd = psd / self.resp

            # Convert to dB
            if self.scale_var == 'linear':
                psd = 10 * np.log10(psd)

            # Interpolate at longer periods
            # mask = np.isnan(psd)
            # psd_interp = np.where(mask, nearest(mask), psd)
            
            # Add to histogram
            self.update_histogram(psd)

        # FCWT method
        elif self.method == 'fcwt':
            iter_obj = self.sliding_window(signal, self.nwin, self.win_shift*self.fs)
            i, Nseg = 0, len(list(iter_obj))
            psd = np.full((Nseg, self.nfreqs), np.nan)

            iter_obj = self.sliding_window(signal, self.nwin, self.win_shift*self.fs)
            for segment in iter_obj:
                _psd, _, _coi = myspec.my_fcwt(self.fcwt_obj, self.morlet, segment, 
                                               self.dt, freqmin=1/self.Tmax, freqmax=1/self.Tmin, 
                                               nptsfreq=self.nfreqs, sampling=self.scale_T, 
                                               smooth=False, mode='psd')
                mask = (self.psd_freqs[:, None] < 1/_coi)
                _psd = np.ma.masked_array(_psd, mask=mask)

                # Average wavelet PSD over time
                psd[i, :] = np.mean(_psd, axis=1)
                i += 1
                
            # Remove instrument response
            psd = psd / self.resp

            # Convert to dB
            if self.scale_var == 'linear':
                psd = 10 * np.log10(psd)
            
            # Add to histogram
            self.update_histogram(psd)

        return psd
        
    def calc_coh(self, signal1, signal2, ns=3, nt=3.0):

        # Welch method
        if self.method == 'welch':
            iter_obj1 = self.sliding_window(signal1, self.nwin, self.win_shift*self.fs)
            i, Nseg = 0, len(list(iter_obj1))
            coh = np.full((Nseg, self.nfreqs), np.nan)

            iter_obj1 = self.sliding_window(signal1, self.nwin, self.win_shift*self.fs)
            iter_obj2 = self.sliding_window(signal2, self.nwin, self.win_shift*self.fs)
            for (seg1, seg2) in zip(iter_obj1, iter_obj2):
                _coh, _ = myspec.my_welch_coh(seg1, seg2, self.dt, self.nperseg, int(0.75*self.nperseg))
                
                # Binning. Use mean for averaging coherence
                for j, per_left, per_right in zip(range(self.nfreqs), self.period_bins[:-1], 
                                                  self.period_bins[1:]):
                    ind = (per_left <= 1/self.psd_freqs) & (1/self.psd_freqs <= per_right)
                    coh[i, j] = np.mean(_coh[ind])
                i += 1

            # Interpolate at longer periods
            # mask = np.isnan(psd)
            # psd_interp = np.where(mask, nearest(mask), psd)
            
            # Add to histogram
            self.update_histogram(coh)

        # FCWT method
        elif self.method == 'fcwt':
            iter_obj1 = self.sliding_window(signal1, self.nwin, self.win_shift*self.fs)
            i, Nseg = 0, len(list(iter_obj1))
            psd = np.full((Nseg, self.nfreqs), np.nan)

            iter_obj1 = self.sliding_window(signal1, self.nwin, self.win_shift*self.fs)
            iter_obj2 = self.sliding_window(signal2, self.nwin, self.win_shift*self.fs)
            for (seg1, seg2) in zip(iter_obj1, iter_obj2):
                _, _, _, _coh, _, _coi = myspec.my_fxwt(self.fcwt_obj, self.morlet, seg1, seg2,
                                                        self.dt, freqmin=1/self.Tmax, freqmax=1/self.Tmin, 
                                                        nptsfreq=self.nfreqs, sampling=self.scale_T, 
                                                        smooth=False, ns=ns, nt=nt)
                mask = (self.psd_freqs[:, None] < 1/_coi)
                _coh = np.ma.masked_array(_coh, mask=mask)

                # Average wavelet PSD over time
                coh[i, :] = np.mean(_coh, axis=1)
                i += 1
            
            # Add to histogram
            self.update_histogram(coh)

        return coh
        
    def update_histogram(self, data):
        """
        Update the histogram based on the given 2D data.

        Parameters:
        - data (array): A 2D numpy array with shape (n_hours, n_periods).
        """
        if self.var_bins is None:
            raise ValueError("Period bins and variable bins must be set before updating the histogram.")

        if self.scale_var not in ["linear", "log"]:
            raise ValueError("Invalid scale. Use 'linear' or 'log'.")

        # Check data consistency with period_bins
        self._check_data(data)

        if self.scale_var == "log":
            data[data <= 0] = np.nan
            data = np.log10(data)

        # Calculate histogram for each period sample
        new_hist = np.apply_along_axis(
            lambda x: np.histogram(x[~np.isnan(x)], bins=self.var_bins)[0], 
            axis=0, arr=data)

        if self.histogram is not None:
            self.histogram += new_hist
        else:
            self.histogram = new_hist

        self.ncount = self.histogram.sum(axis=0).max()

    def reduce_var_bins(self, factor=2, inplace=False):
        """
        Reduce the number of bins by a factor.

        Parameters:
        - factor (int): Reduce the number of bins by this factor.
        """

        # Create new variable bins
        new_var_bins = self.var_bins[::factor]
        if new_var_bins[-1] < self.var_bins[-1]:
            if self.scale_var == 'log':
                new_var_bins = np.append(new_var_bins, new_var_bins[-1]**2 / new_var_bins[-2])

            elif self.scale_var == 'linear':
                new_var_bins = np.append(new_var_bins, new_var_bins[-1]*2 - new_var_bins[-2])
        new_var_bins_centers = self.convert_to_bin_centers(new_var_bins, self.scale_var)

        # Sum the counts in adjacent bins
        new_histogram = np.add.reduceat(self.histogram, np.arange(0, self.histogram.shape[0], factor), axis=0)

        if inplace:
            self.var_bins, self.var_bins_centers = new_var_bins, new_var_bins_centers
            self.histogram = new_histogram
        else:
            return new_var_bins, new_histogram

    def _check_data(self, data):
        """
        Check if the shape of the input data is consistent with period_bins.

        Parameters:
        - data (array): 2D numpy array with shape (n_hours, n_periods).
        """
        if data.shape[1] != len(self.period_bins) - 1:
            raise ValueError("The number of period samples for the input data does not match period_bins.")

    def __str__(self):
        """
        String representation of the myPDF object.
        """
        if self.histogram is None:
            return "Empty histogram."
        
        _period_bins = [f"{val:.2f}" for val in self.period_bins]
        _var_bins = [f"{val:.2f}" for val in self.var_bins]

        var_bin_labels = [f"{_var_bins[i]}-{_var_bins[i+1]}" \
                          for i in range(len(self.var_bins)-1)]
        period_labels = [f"{_period_bins[i]}-{_period_bins[i+1]} s" \
                         for i in range(len(self.period_bins)-1)]

        df = pd.DataFrame(self.histogram.T, columns=var_bin_labels, index=period_labels)
        return f"Scale_T: {self.scale_T}\nScale_var: {self.scale_var}\nPeriod Bins: {_period_bins}\nVariable Bins: {_var_bins}\nHistogram:\n{df}"
    
    def get_percentile(self, percentile=50):
        """
        Calculate the percentile at each period sample from the histogram.

        Parameters:
        - percentile (float, optional): The percentile to calculate. Default is 50.

        Returns:
        - pc_hist (float): The percentile to obtain.
        """
        if self.histogram is None:
            raise ValueError("Histogram is empty. Please update the histogram first.")

        # Calculate the cumulative sum
        pdf_hist = self.histogram / np.sum(self.histogram, axis=0, keepdims=True)
        cdf_hist = np.cumsum(pdf_hist, axis=0)

        # Find the index of the first value at each period sample
        index_pc = np.argmax(cdf_hist >= percentile/100, axis=0)

        # Get the corresponding bin values
        pc_hist = np.array([self.var_bins_centers[idx] for idx in index_pc])

        # Remove empty bins
        mask = np.isnan(pdf_hist[0])
        pc_hist = np.ma.masked_array(pc_hist, mask=mask)

        return pc_hist
    
    def get_mode(self):
        """
        Get the mode values of the histogram.

        Returns:
        - mode_hist (array): An array containing the mode values at each period sample.
        """
        if self.histogram is None:
            raise ValueError("Histogram is empty. Please update the histogram first.")

        # Get the index of the maximum count (mode) at each period sample
        index_mode = np.argmax(self.histogram, axis=0)

        # Get the corresponding bin values
        mode_hist = np.array([self.var_bins_centers[idx] for idx in index_mode])

        # Remove empty bins
        mask = np.isnan(self.histogram[0])
        mode_hist = np.ma.masked_array(mode_hist, mask=mask)

        return mode_hist

    def plot_histogram(self, ylabel=None, vmin=None, vmax=None, cmap=None, 
                       ticks_loc=None, show_noise_models=True):
        """
        Plot the histogram with probability using pcolormesh with flat shading.

        Parameters:
        - ylabel (str, optional): The label for the y-axis. Default is None.
        - vmax (float, optional): The maximum value for the color scale. Default is None.
        - cmap (str, optional): The colormap for the plot. Default is None.
        """
        if self.histogram is None:
            raise ValueError("Histogram is empty. Please update the histogram first.")
            
        # Convert to probability
        pdf_hist = self.histogram / np.sum(self.histogram, axis=0, keepdims=True) \
            / np.diff(self.var_bins)[:, None]

        # Remove empty bins
        mask = ~np.isnan(pdf_hist[0])
        xdata = self.period_bins[np.concatenate([[True], mask])]
        
        # Color range for probability
        if vmin is None:
            vmin = 0
        if vmax is None:
            vmax = np.nanmax(pdf_hist)
            
        # Colormap
        if cmap is None:
            cmap = 'jet'

        # Plot the histogram using pcolormesh with flat shading
        fig = plt.figure(figsize=(6,3), dpi=200)
        ax = fig.add_subplot()
        obj = ax.pcolormesh(xdata, self.var_bins, pdf_hist[:, mask], 
                            vmin=vmin, vmax=vmax, shading='flat', cmap=cmap)

        # Set axis properties
        if ylabel is None:
            ylabel = r'PSD dB (10$\times$log$_{10}$(m$^2$/s$^4$/Hz))'
        ax.set_xscale(self.scale_T)
        ax.set_yscale(self.scale_var)
        ax.set_xlabel('Period (s)')
        ax.set_ylabel(ylabel)
        
        # Set ticks
        if ticks_loc is None:
            ticks_loc = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 
                         100, 200, 500]
        ax.set_xticks(ticks_loc, ticks_loc)
        ax.grid()

        # Show colorbar
        fig.colorbar(obj, label='Probability Density', ax=ax)

        # Plot noise models
        if show_noise_models:
            models = (get_nhnm(), get_nlnm())
            for periods, noise_model in models:
                ax.plot(periods, noise_model, '1.0', linewidth=1.0, zorder=10)

        # Restrict to relevant domain
        ax.set_xlim(xdata[0], xdata[-1])
        ax.set_ylim(self.var_bins[0], self.var_bins[-1])

        # Print number of PSD curves
        print(f'Total {self.ncount} PSD curves.')

        return fig, ax

    def save_npz(self, filename):
        """
        Saves the PPSD as a compressed numpy binary (npz format).
        The resulting file can be restored using `my_ppsd.load_npz(filename)`.

        Parameters:
        - filename (str): Name of numpy .npz output file
        """

        out = {}
        for key in (self.NPZ_STORE_KEYS_BASIC + self.NPZ_STORE_KEYS_PSD):
            value = getattr(self, key)
            if key in self.NPZ_STORE_KEYS_FIX:
                try:
                    value = self.NPZ_SIMPLE_TYPE_MAP.get(value, value)
                except:
                    pass
            out[key] = value
        np.savez_compressed(filename, **out)

    @staticmethod
    def load_npz(filename):
        """
        Load previously computed PSD results from a
        compressed numpy binary in npz format.

        Parameters:
        - filename (str): Name of numpy .npz file with stored PPSD data
        """
        def _load(data):

            # Load stored data into myPDF class
            mypdf = myPDF()

            for key in (myPDF.NPZ_STORE_KEYS_BASIC + myPDF.NPZ_STORE_KEYS_PSD):
                try:
                    data_ = data[key]
                except ValueError:
                    msg = ("The allow_pickle parameter of "
                           "myPDF.load_npz is always set to True")
                    raise ValueError(msg)
                
                if key in myPDF.NPZ_STORE_KEYS_STR_TYPES:
                    data_ = str(data_)
                elif key in (myPDF.NPZ_STORE_KEYS_FIX):
                    try:
                        data_ = myPDF.NPZ_SIMPLE_TYPE_MAP_R.get(data_, data_)
                    except:
                        pass
                setattr(mypdf, key, data_)

            return mypdf
        
        with np.load(filename, allow_pickle=True) as data:
            mypdf = _load(data)

            # Calculate relevant parameters
            mypdf.set_data_stats(mypdf.channel, mypdf.window, mypdf.resp_comp)
            mypdf.period_bins_centers = mypdf.convert_to_bin_centers(mypdf.period_bins, mypdf.scale_T)
            mypdf.var_bins_centers = mypdf.convert_to_bin_centers(mypdf.var_bins, mypdf.scale_var)
            mypdf.ncount = mypdf.histogram[:,2].sum()

            # Initialize FFT analyzer
            mypdf.init_fft_analyzer(mypdf.method, mypdf.win_shift, mypdf.nperseg, mypdf.w0)
                
            return mypdf
        
    def trim_trace(self, trace):
            
        L = 86400 * self.fs         # Exact one-day data samples
        npts = trace.stats.npts     # Number of data points
        X = 3600                    # 1-hour in seconds
        ts = trace.stats.starttime  # Trace start time

        # Qualified trace
        if np.abs(npts - L) <= 10 * self.fs:
            return 0

        # Find the nearest hour
        if ts.timestamp % X <= 60:
            return ts.hour
        else:
            ts_ = UTCDateTime(np.ceil(ts.timestamp / X) * X)
            trace.trim(starttime=ts_)
            return ts.hour+1


### Noise models (Peterson, 1993) ###
NOISE_MODEL_FILE = '/Users/qingji/Documents/Datasets/noise_models/noise_models.npz'

def get_nlnm():
    data = np.load(NOISE_MODEL_FILE)
    periods = data['model_periods']
    nlnm = data['low_noise']
    return (periods, nlnm)

def get_nhnm():
    data = np.load(NOISE_MODEL_FILE)
    periods = data['model_periods']
    nhnm = data['high_noise']
    return (periods, nhnm)