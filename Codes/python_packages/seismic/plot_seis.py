#!/usr/bin/env python3

"""
Plotting functions for seismic analysis

Author: Qing Ji
"""

# Load python packages
import fnmatch
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Plot one trace
def plot_one_trace(trace_in, filt=None, filt_param={'corners':4, 'zerophase':True}, 
                   taper=None, ax=None, xaxis_date=True, title=None, **kwargs):

    # ObsPy Trace object
    trace = trace_in.copy()

    # Filter
    if filt is not None:
        trace.detrend()
        if taper is not None:
            trace.taper(taper)
        if filt[0] is None:
            trace.filter(type='lowpass', freq=filt[1], **filt_param)
        elif filt[1] is None:
            trace.filter(type='highpass', freq=filt[0], **filt_param)
        else:
            trace.filter(type='bandpass', freqmin=filt[0], freqmax=filt[1], **filt_param)

    # Create new Axes object
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 2))
    
    # Time axis & Plot data
    if xaxis_date:
        ax.plot(trace.times("matplotlib"), trace.data, "-", **kwargs)
        ax.xaxis_date()
        ax.set_xlim(trace.times("matplotlib")[0], trace.times("matplotlib")[-1])
    else:
        ax.plot(trace.times(), trace.data, "-", **kwargs)
        ax.set_xlabel('Time (s)')
        ax.set_xlim(trace.times()[0], trace.times()[-1])

    # Title
    if title is None:
        ax.set_title('%s' % ('.'.join((trace.id).split('.')[:2])))
    else:
        ax.set_title(title)

    # Y-axis label
    if hasattr(trace.stats, 'unit'):
        ax.set_ylabel(trace.stats.unit)
    
    if 'fig' in locals():
        return [fig, ax]
    else:
        return ax
    

# Plot two records in the same panel
def plot_two_traces(trace1_in, trace2_in, filt=None, **kwargs):

    trace1 = trace1_in.copy()
    trace2 = trace2_in.copy()
    # Band-pass filt
    if filt is not None:
        trace1.filter(type='bandpass', freqmin=filt[0], freqmax=filt[1])
        trace2.filter(type='bandpass', freqmin=filt[0], freqmax=filt[1])

    fig, ax1 = plt.subplots(figsize=(20, 5))
    ax1.plot(trace1.times("matplotlib"), trace1.data,
             "k-", **kwargs)
    ax1.xaxis_date()
    ax1.set_title('%s' % ('.'.join((trace1.id).split('.')[:2])))
    low, high = ax1.get_ylim()
    bound = max(abs(low), abs(high))
    ax1.set_ylim(-bound, bound)

    ax_dist = ax1.twinx()
    ax_dist.plot(trace2.times("matplotlib"), trace2.data, "r-",
                 alpha=0.9, **kwargs)
    ax_dist.tick_params(axis='y', colors='r')
    low, high = ax_dist.get_ylim()
    bound = max(abs(low), abs(high))
    ax_dist.set_ylim(-bound, bound)

    if fnmatch.fnmatch(trace1.stats.channel, '*D[FH]'):
        ax1.set_ylabel('Pressure (Pa)')
    else:
        ax1.set_ylabel('Displacement ($\mu$m)')

    if fnmatch.fnmatch(trace2.stats.channel, '*D[FH]'):
        ax_dist.set_ylabel('Pressure (Pa)', color='r')
    else:
        ax_dist.set_ylabel('Displacement ($\mu$m)', color='r')

    return fig, [ax1, ax_dist]


def plot_map_data(extent, lon, lat, proj=None, line=False, ax=None, **kwargs):

    # Map projection
    if proj is None:
        proj = ccrs.Mercator()
    
    # Data projection (lon/lat input data)
    proj_data = ccrs.PlateCarree()

    # Figure setup
    flag = False
    if ax is None:
        fig = plt.figure()
        ax = plt.axes(projection=proj)
        ax.set_extent(extent, crs=proj_data)
        ax.add_feature(cfeature.STATES, edgecolor='gray')
        ax.add_feature(cfeature.COASTLINE, edgecolor='gray')
        flag = True

    # Plot data points
    if line:
        ax.plot(lon, lat, transform=proj_data, **kwargs)
    else:
        ax.scatter(lon, lat, transform=proj_data, **kwargs)

    # Grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--')
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    gl.top_labels = False
    gl.right_labels = False

    if flag:
        return fig, ax