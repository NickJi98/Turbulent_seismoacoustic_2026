#!/usr/bin/env python3

"""
Plotting functions for Fourier analysis

Author: Qing Ji
"""

# Load python packages
import numpy as np
import fnmatch

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from cmcrameri import cm


# Plot wavelet transform results
# Defualt figure properties
cb_pad_wide, cb_pad_short = 0.08, 0.04
cb_pad_wide, cb_pad_short = 0.11, 0.04
line_width = 1.0


# Plot wavelet PSD
def plot_wt_psd(timestamp, freqs, cfs, coi=None, channel='LHZ', decimate=50,
                cmap=cm.vik, freq_scale='log', value_range=[-4,4], title='Title',
                tc_time=None, tc_dist=None, dist_range=[0,1e3], dist_color='white',
                cb_label=None, figsize=(9,1.8), ax=None, xaxis_date=True):
    
    # Colorbar location
    if tc_time is None:
        cb_pad = cb_pad_short
    else:
        cb_pad = cb_pad_wide
   
    # Figure & Axes objects
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Wavelet PSD
    obj = ax.pcolormesh(timestamp[::decimate], freqs, np.log10(cfs[:, ::decimate]), 
                        cmap=cmap, edgecolors='none',
                        vmin=value_range[0], vmax=value_range[1])
    if channel is None:
        label = cb_label
    else:
        if fnmatch.fnmatch(channel, '*D[FH]'):
            label = 'Log Pressure PSD (Pa$^2$/Hz)'
            title = 'Wavelet Pressure PSD, %s' %(title)
        elif fnmatch.fnmatch(channel, '*HZ'):
            label = 'Log Disp. PSD ($\mu$m$^2$/Hz)'
            title = 'Wavelet Vertical Disp. PSD, %s' %(title)
        elif fnmatch.fnmatch(channel, '*HN'):
            label = 'Log Disp. PSD ($\mu$m$^2$/Hz)'
            title = 'Wavelet N-S Disp. PSD, %s' %(title)
        elif fnmatch.fnmatch(channel, '*HE'):
            label = 'Log Disp. PSD ($\mu$m$^2$/Hz)'
            title = 'Wavelet E-W Disp. PSD, %s' %(title)
    cb = fig.colorbar(obj, ax=ax, pad=cb_pad, label=label)
    cb.outline.set_linewidth(0.5)
   
    # Station distance from hurricane
    if tc_time is not None:
        ax_dist = ax.twinx()
        ax_dist.plot(tc_time, tc_dist, color=dist_color, linewidth=line_width)
        ax_dist.tick_params(axis='y', colors='r')
        ax_dist.set_ylabel('Distance (km)', color='r')
        ax_dist.set_ylim(dist_range[0], dist_range[1])
   
    # Cone of influence (COI)
    if coi is not None:
        ax.plot(timestamp[::decimate], 1/coi[::decimate], 'w--', linewidth=line_width)

    if xaxis_date:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(
            mdates.AutoDateLocator(minticks=6, maxticks=10)))
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
    ax.set_xlim(timestamp[0], timestamp[-1])
    ax.set_ylabel('Frequency (Hz)')
    ax.set_yscale(freq_scale)
    ax.set_ylim(freqs[-1], freqs[0])
    ax.set_title(title)

    if tc_time is not None:
        return [fig, ax, ax_dist]
    else:
        return [fig, ax]


# Plot wavelet CSD amplitude
def plot_wt_csd_amp(timestamp, freqs, WX, coi=None, decimate=50,
                    cmap=cm.vik, freq_scale='log', title='Title',
                    tc_time=None, tc_dist=None, dist_range=[0,1e3], 
                    dist_color='white', xaxis_date=True, figsize=(9,1.8), ax=None):

    # Colorbar location
    if tc_time is None:
        cb_pad = cb_pad_short
    else:
        cb_pad = cb_pad_wide

    # Figure & Axes objects
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # CSD amplitude
    obj = ax.pcolormesh(timestamp[::decimate], freqs, np.log10(np.abs(WX[:, ::decimate])), 
                        cmap=cmap, edgecolors='none')
    cb = fig.colorbar(obj, ax=ax, pad=cb_pad, label='Log CSD Amplitude')
    cb.outline.set_linewidth(0.5)

    # Station distance from hurricane
    if tc_time is not None:
        ax_dist = ax.twinx()
        ax_dist.plot(tc_time, tc_dist, color=dist_color, linewidth=line_width)
        ax_dist.tick_params(axis='y', colors='r')
        ax_dist.set_ylabel('Distance (km)', color='r')
        ax_dist.set_ylim(dist_range[0], dist_range[1])
        
    # Cone of influence (COI)
    if coi is not None:
        ax.plot(timestamp[::decimate], 1/coi[::decimate], 'w--', linewidth=line_width)
        
    if xaxis_date:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(
            mdates.AutoDateLocator(minticks=6, maxticks=10)))
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
    ax.set_xlim(timestamp[0], timestamp[-1])
    ax.set_ylabel('Frequency (Hz)')
    ax.set_yscale(freq_scale)
    ax.set_ylim(freqs[-1], freqs[0])
    ax.set_title('Wavelet CSD Amplitude, %s' %title)

    if tc_time is not None:
        return [fig, ax, ax_dist]
    else:
        return [fig, ax]


# Plot wavelet CSD phase
def plot_wt_csd_phase(timestamp, freqs, WX, coi=None, decimate=50,
                      cmap=cm.vik, freq_scale='log', title='Title',
                      tc_time=None, tc_dist=None, dist_range=[0,1e3], 
                      dist_color='white', xaxis_date=True, figsize=(9,1.8), ax=None):
    
    # Colorbar location
    if tc_time is None:
        cb_pad = cb_pad_short
    else:
        cb_pad = cb_pad_wide

    # Figure & Axes objects
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # CSD phase
    obj = ax.pcolormesh(timestamp[::decimate], freqs, np.rad2deg(np.angle(WX[:, ::decimate])), 
                        cmap=cmap, edgecolors='none', vmin=-180, vmax=180)
    cb = fig.colorbar(obj, ax=ax, pad=cb_pad, label='Phase Diff. (deg)', ticks=[-180, -90, 0, 90, 180])
    cb.outline.set_linewidth(0.5)

    # Station distance from hurricane
    if tc_time is not None:
        ax_dist = ax.twinx()
        ax_dist.plot(tc_time, tc_dist, color=dist_color, linewidth=line_width)
        ax_dist.tick_params(axis='y', colors='r')
        ax_dist.set_ylabel('Distance (km)', color='r')
        ax_dist.set_ylim(dist_range[0], dist_range[1])

    # Cone of influence (COI)
    if coi is not None:
        ax.plot(timestamp[::decimate], 1/coi[::decimate], 'w--', linewidth=line_width)

    if xaxis_date:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(
            mdates.AutoDateLocator(minticks=6, maxticks=10)))
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
    ax.set_xlim(timestamp[0], timestamp[-1])
    ax.set_ylabel('Frequency (Hz)')
    ax.set_yscale(freq_scale)
    ax.set_ylim(freqs[-1], freqs[0])
    ax.set_title('Wavelet CSD Phase, %s' %title)

    if tc_time is not None:
        return [fig, ax, ax_dist]
    else:
        return [fig, ax]


# Plot coherence of two traces
def plot_wt_coh(timestamp, freqs, Wcoh, coi=None, decimate=50,
                cmap=cm.vik, freq_scale='log', value_range=[0,1], title='Title',
                tc_time=None, tc_dist=None, dist_range=[0,1e3], 
                dist_color='white', xaxis_date=True, levels=None,
                figsize=(9,1.8), ax=None):
    
    # Colorbar location
    if tc_time is None:
        cb_pad = cb_pad_short
    else:
        cb_pad = cb_pad_wide

    # Figure & Axes objects
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Coherence
    obj = ax.pcolormesh(timestamp[::decimate], freqs, Wcoh[:, ::decimate],
                        cmap=cmap, edgecolors='none',
                        vmin=value_range[0], vmax=value_range[1])

    # Contours
    if levels is not None:
        ax.contour(timestamp[::decimate], freqs, Wcoh[:, ::decimate], levels)
    
    cb = fig.colorbar(obj, ax=ax, pad=cb_pad, label='Coherence')
    cb.outline.set_linewidth(0.5)

    # Station distance from hurricane
    if tc_time is not None:
        ax_dist = ax.twinx()
        ax_dist.plot(tc_time, tc_dist, color=dist_color, linewidth=line_width)
        ax_dist.tick_params(axis='y', colors='r')
        ax_dist.set_ylabel('Distance (km)', color='r')
        ax_dist.set_ylim(dist_range[0], dist_range[1])

    # Cone of influence (COI)
    if coi is not None:
        ax.plot(timestamp[::decimate], 1/coi[::decimate], 'w--', linewidth=line_width)

    if xaxis_date:
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(
            mdates.AutoDateLocator(minticks=6, maxticks=10)))
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))
    ax.set_xlim(timestamp[0], timestamp[-1])
    ax.set_ylabel('Frequency (Hz)')
    ax.set_yscale(freq_scale)
    ax.set_ylim(freqs[-1], freqs[0])
    ax.set_title('Coherence, %s' %(title))

    if tc_time is not None:
        return [fig, ax, ax_dist]
    else:
        return [fig, ax]


# Plot PSD time evolution at particular frequencies
def plot_time_evolution(timestamp, freqs, cfs, freq_samples, channel,
                        title='Title', legend_unit='period'):

    fig, ax = plt.subplots(figsize=(20, 5))
    if fnmatch.fnmatch(channel, '*D[FH]'):
        label = 'Pressure PSD [Pa$^2$/Hz]'
        title = 'Wavelet Pressure PSD, %s' %title
    else:
        label = 'Disp. PSD [($\mu$m)$^2$/Hz]'
        title = 'Wavelet Vertical Disp. PSD, %s' %title

    for freq_point in freq_samples:
        freq_ind = (np.abs(freqs - freq_point)).argmin()
        freq_trace = cfs[freq_ind, :]
        if legend_unit == 'period':
            ax.plot(timestamp, freq_trace, label='%d s' %(1/freq_point))
        else:
            ax.plot(timestamp, freq_trace, label='%.2f Hz' %freq_point)

    ax.set_xlim(timestamp[0], timestamp[-1])
    ax.set_ylabel(label)
    ax.set_title(title)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=5)
    ax.grid()
    fig.show()

    return [fig, ax]


# Plot PSD snapshot
def plot_psd_snapshot(psd_db, comp='prs', snapshot_time=None, 
                      az_color=False, show_cb=True):
    
    if show_cb:
        fig, ax = plt.subplots(figsize=(10, 10))
    else:
        fig, ax = plt.subplots(figsize=(10, 7))
        
    # Inter-quartile range over 1-hr
    ax.errorbar(psd_db['dist'], psd_db[f'{comp}_med'], 
                yerr=[psd_db[f'{comp}_quar1'], psd_db[f'{comp}_quar3']], 
                fmt='none', c='k', alpha=0.3)
    
    # Median PSD level
    if az_color:
        obj = ax.scatter(psd_db['dist'], psd_db[f'{comp}_med'], 
                         c=psd_db['azimuth'].tolist(), 
                         cmap='hsv', vmin=-180, vmax=180)
    else:
        ax.scatter(psd_db['dist'], psd_db[f'{comp}_med'], c='k')
        
    if comp == 'prs':
        ax.set_ylabel('Pressure PSD [Pa$^2$/Hz]')
        ax.set_ylim([1e-3, 1e4])
        
    elif comp == 'uz':
        ax.set_ylabel('Displacement PSD [($\mu$m)$^2$/Hz]')
        ax.set_ylim([1e-4, 1e3])
        
    ax.set_xlabel('Distance from Hurricane Center [km]')
    ax.set_yscale('log')
    ax.set_xlim([0, 1e3])
    ax.grid()
    
    # Show colorbar for azimuthal range
    if az_color and show_cb:
        fig.colorbar(obj, label='Azimuth [deg]', orientation='horizontal')
    
    if snapshot_time is not None:
        ax.set_title(snapshot_time.strftime('%m-%d %H:%M'))
        
    return fig, ax