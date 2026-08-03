"""
Common paths and packages
Load this file before running the notebooks
"""

### Directories
import os
main_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Parent directory
data_dir = os.path.join(main_dir, 'Data')
fig_dir = os.path.join(main_dir, 'Figures/Panels')
sup_dir = os.path.join(main_dir, 'Figures_Sup/Panels')


### Python packages
import sys, glob, copy, time
import numpy as np
import pandas as pd
import netCDF4 as nc

# scipy
from scipy.io import loadmat
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from scipy.integrate import cumulative_trapezoid
from scipy.stats import circmean, circstd
from scipy.optimize import leastsq, curve_fit
from scipy.interpolate import griddata
from scipy import fft

# Plotting
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# Apply the shared plotting style
mpl.rc_file(os.path.join(main_dir, 'Codes/python_packages/matplotlibrc'))
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cv2
from pyproj import Transformer, Geod

# xarray dataset
import xarray as xr
from xgrads import CtlDescriptor, open_CtlDataset

# obspy
from obspy import Stream, UTCDateTime, read
from obspy.clients.fdsn import Client
from datetime import datetime, timedelta

# Temporary fix for using fcwt on Mac (likely libomp.dylib issue)
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


### Python functions under {main_dir}/Codes
codes_dir = os.path.join(main_dir, "Codes/python_packages/")
if codes_dir not in sys.path:
    sys.path.append(codes_dir)

from seismic.seis_download import download_trace, check_npts
import Fourier.spectrum as fspec
import Fourier.plot_spec as pspec
import Fourier.MISO as miso
import TC_analysis.atm_dataset as atm_ds
import TC_analysis.tc_process as tc_func
import TC_analysis.cm1_analysis as cm1_func


### Some station location
rs_sta_loc = np.array([-89.8250, 30.3369])          # Radiosonde ID 72233
tower_loc = np.array([-90.693993, 29.648682])       # Wind tower T3, Isaac 2012


### Functions: Read data and info ###

def get_sta_loc(sta_char, sta_info):
    
    _mask = (sta_info['net'] == sta_char[0]) & (sta_info['stnm'] == sta_char[1])
    sta_loc = np.squeeze(sta_info[_mask][['lon', 'lat']].to_numpy())
    
    return sta_loc

def get_asos_sta_loc(asos_char, asos_info):
    
    _mask = (asos_info['stid'] == asos_char)
    sta_loc = np.squeeze(asos_info[_mask][['lon', 'lat']].to_numpy())
    
    return sta_loc

def read_asos_metadata(csv_files, ts=None, te=None):
    
    # Read each network metadata
    dfs = []
    for file in csv_files:
        df = pd.read_csv(file)
        dfs.append(df)

    # Combine all networks
    asos_info = pd.concat(dfs, ignore_index=True)
    asos_info['begints'] = pd.to_datetime(asos_info['begints'])
    asos_info['endts'] = pd.to_datetime(asos_info['endts'])
    
    # Select based on time range
    if (ts is None) and (te is None):
        asos_filt_info = asos_info
    elif te is None:
        asos_filt_info = asos_info[asos_info['begints'] <= ts]
    else:
        asos_filt_info = asos_info[(asos_info['begints'] <= ts) & ((asos_info['endts'] >= te) | asos_info['endts'].isnull())]

    return asos_filt_info

def read_tower_data(file_list, clean_data=True, clean_thresh=0.01, print_info=True):

    # Single file input
    if isinstance(file_list, str):
        file_list = [file_list]
	
    # Timestamp of start time & end time
    mat_data = loadmat(file_list[0])
    time_start = datetime.strptime(mat_data['Time_Start_10m'][0], '%d-%b-%Y %H:%M:%S')
    mat_data = loadmat(file_list[-1])
    time_end = datetime.strptime(mat_data['Time_End_10m'][0], '%d-%b-%Y %H:%M:%S')

    if print_info:
        print("Start time:", time_start)
        print("End time:", time_end)
        print()

    # Basic time axis
    dt = 1 / 32
    t_ = np.arange(28800) * dt

    datasets = []
    for i, file in enumerate(file_list):

        # Extract variables for all heights
        heights = ['Height_5m', 'Height_10m', 'Height_15m',
                   'Height_75m', 'Height_125m']
        data_dict = {}

        # Time axis
        t = t_ + i * 900

        mat_data = loadmat(file)
        for height in heights:
            ws = mat_data['USGill_Speed'][0,0][height].flatten()        # Wind speed
            wd = mat_data['USGill_Direction'][0,0][height].flatten()    # Wind direction
            w = mat_data['USGill_Vertical'][0,0][height].flatten()      # Vertical wind
            T = mat_data['USGill_Temperature'][0,0][height].flatten()   # Temperature
            
            # Convert wind direction to u/v components
            wd_rad = np.deg2rad(wd)
            u1 = -ws * np.sin(wd_rad)   # West-to-East
            u2 = -ws * np.cos(wd_rad)   # South-to-North
            
            # Store in dictionary
            data_dict[height] = xr.Dataset(
                data_vars={
                    'u1': (['time'], u1),
                    'u2': (['time'], u2),
                    'w': (['time'], w),
                    'T': (['time'], T),
                    'ws': (['time'], ws),
                    'wd': (['time'], wd)
                },
                coords={'time': t}
            )
        
        # Combine into a single Dataset with 'height' as a dimension
        ds = xr.concat(data_dict.values(), dim='height')
        ds['height'] = [int(h.split('_')[-1].replace('m','')) for h in heights]
        ds['height'] = ds['height'].where(~ds['height'].isin([75, 125]), ds['height'] / 10)
        datasets.append(ds)

    # Concatenate all 15-minute files along the time dimension
    raw_ds = xr.concat(datasets, dim='time').sortby('time')
    Nt = raw_ds.dims['time']

    # Add time_start and time_end as global attributes
    raw_ds.attrs['time_start'] = time_start.strftime("%Y-%m-%dT%H:%M:%S")
    raw_ds.attrs['time_end'] = time_end.strftime("%Y-%m-%dT%H:%M:%S")


    ### Clean NaN values & Interpolation ###
    if clean_data:
        if print_info:
            print("NaN counts by height and variable:")
        nan_ds = raw_ds.isnull()

        for height in nan_ds.height.values:
            # Get NaN counts for each height
            _tmp = nan_ds.sel(height=height).sum()
            _str = ', '.join([f"{var}: {_tmp[var].item()}" 
                                    for var in _tmp.data_vars])
            if print_info:
                print(f"{height}m: {_str}")
    
        # Remove height levels with too many NaNs
        keep_heights = []
        if print_info:
            print()

        for height in nan_ds.height.values:
            _tmp = nan_ds.sel(height=height).sum()
            nan_frac = np.mean([_tmp[var].item() for var in _tmp.data_vars]) / Nt
            if nan_frac <= clean_thresh:
                keep_heights.append(height)
            else:
                if print_info:
                    print(f"Excluding {height}m: {nan_frac*100:.2f}% NaNs")

        wind_ds = raw_ds.sel(height=keep_heights)
        wind_ds = wind_ds.interpolate_na(dim='time', method='linear', fill_value="extrapolate")
        if print_info:
            print("Height levels retained:", keep_heights)
        return wind_ds
    
    else:
        return raw_ds
    
def read_dropsonde_metadata():

    # Dropsonde files
    nc_filelist = glob.glob(os.path.join(data_dir, 'Isaac/dropsonde/', '20120828H2/', '*.nc')) \
            + glob.glob(os.path.join(data_dir, 'Isaac/dropsonde/', '20120829A1/', '*.nc'))

    # Create Dataframe
    time_list, lon_list, lat_list = [], [], []

    # Obtain time of the dropsonde
    def get_sonde_reftime(nc_db):
        time_str = nc_db['launch_time'].units
        _, date_time_str = time_str.split('since')
        return UTCDateTime(date_time_str[:-4])

    for nc_file in nc_filelist:
        nc_db = nc.Dataset(nc_file, 'r')
        ref_time = get_sonde_reftime(nc_db)
        ref_lon = nc_db['reference_lon'][:][0]
        ref_lat = nc_db['reference_lat'][:][0]
        
        if np.ma.is_masked(ref_lon):
            ref_lon = nc_db['lon'][:].compressed()[0]
        if np.ma.is_masked(ref_lat):
            ref_lat = nc_db['lat'][:].compressed()[0]
        
        time_list.append(ref_time)
        lon_list.append(ref_lon)
        lat_list.append(ref_lat)
        
    nc_df = pd.DataFrame({'ref_time': time_list, 'ref_lon': lon_list, 'ref_lat': lat_list, 
                        'filename': nc_filelist})
    return nc_df


### Function: Structure function ###

def calc_Dst(u, dt, window_length, step=None, max_lag=None):
    """
    Compute second-order structure function D(tau) using overlapping windows.

    Parameters
    ----------
    u : ndarray
        Input time series.
    dt : float
        Sampling interval [s].
    window_length : float
        Window length [s].
    step : float or None
        Step size between windows [s]. Smaller than window_length for overlap.
        If None, defaults to half the window_length (50% overlap).
    max_lag : float or None
        Maximum lag to compute [s]. If None, defaults to half window length.

    Returns
    -------
    tau : ndarray
        Time lags [s].
    Dst : ndarray
        Structure function values (averaged across windows).
    """

    nwin = int(window_length / dt)
    npts = len(u)

    if step is None:
        step = window_length
    nstep = int(step / dt)

    if max_lag is None:
        max_lag = window_length / 2
    nlag = int(max_lag / dt)

    Dst_all = []

    # Slide with overlap
    for start in range(0, npts - nwin + 1, nstep):
        segment = u[start:start+nwin]
        diffs = np.zeros(nlag)

        for lag in range(1, nlag):
            du = segment[lag:] - segment[:-lag]
            diffs[lag] = np.mean(du**2)

        Dst_all.append(diffs)

    Dst_all = np.array(Dst_all)
    Dst = np.nanmean(Dst_all, axis=0)

    tau = np.arange(nlag) * dt
    return Dst, tau


### Functions: Plotting ###

def make_basemap(tc_track_tf=None, tc_info=None, legend=False, date_marker=False, t0=None, f_loc=None, M=None):
    
    # Figure size
    FIGSIZE = 13.333
    img_range = (0, 26.6667, 0, 14.8333)
    
    # Basemap
    map_img = plt.imread(os.path.join(data_dir, 'Isaac_HQ.png'))
    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE * map_img.shape[0] / map_img.shape[1]), dpi=300)
    ax.imshow(map_img, zorder=0, extent=img_range)
    
    # Plot TC track
    if tc_track_tf is not None and tc_info is not None:
        marker_size = np.sqrt(100)
        track_width = 3
        obj_tc_track, = ax.plot(tc_track_tf[:, 0], tc_track_tf[:, 1], 'b-',
                                label='Hurricane Track', zorder=20, lw=track_width)
        obj_dot_1, = ax.plot(tc_track_tf[:, 0][tc_info['scale'] >= 1],
                             tc_track_tf[:, 1][tc_info['scale'] >= 1],
                             ms=marker_size, c='yellow', marker='o', mec='None', ls='None',
                             label='Category 1 Hurricane', zorder=21)
        obj_dot_2, = ax.plot(tc_track_tf[:, 0][tc_info['scale'] == 0],
                             tc_track_tf[:, 1][tc_info['scale'] == 0],
                             ms=marker_size, c='#4DFFFF', marker='o', mec='None', ls='None',
                             label='Tropical Storm', zorder=21)
        obj_dot_3, = ax.plot(tc_track_tf[:, 0][tc_info['scale'] == -1],
                             tc_track_tf[:, 1][tc_info['scale'] == -1],
                             ms=marker_size, c='#6EC1EA', marker='o', mec='None', ls='None',
                             label='Tropical Depression', zorder=21)
    
    # Plot date marker
    if date_marker and t0 is not None:
        for tc_time in [UTCDateTime('20120828-00'), UTCDateTime('20120829-00'), UTCDateTime('20120830-00')]:
            tc_lon = f_loc[0](tc_time-UTCDateTime(t0))
            tc_lat = f_loc[1](tc_time-UTCDateTime(t0))
            tc_loc_tf = np.squeeze(cv2.perspectiveTransform(np.float32([tc_lon, tc_lat])[None, None], M))
            obj_date, = ax.plot(tc_loc_tf[0], tc_loc_tf[1], ms=marker_size, c='r', marker='o', mec='None', ls='None', 
                                zorder=23, label='Date Marker')

        for tc_time in [UTCDateTime('20120828-00'), UTCDateTime('20120829-00'), UTCDateTime('20120830-00')]:
            tc_lon = f_loc[0](tc_time-UTCDateTime(t0))
            tc_lat = f_loc[1](tc_time-UTCDateTime(t0))
            tc_loc_tf = np.squeeze(cv2.perspectiveTransform(np.float32([tc_lon, tc_lat])[None, None], M))
            obj_text = ax.text(tc_loc_tf[0]+0.1, tc_loc_tf[1]+0.1, r'%02d/%02d' %(tc_time.month, tc_time.day), 
                               color='red', horizontalalignment='left', verticalalignment='bottom', 
                               fontsize=20, fontweight='bold', zorder=22)
            obj_text.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='None'))
    
    ax.set_xlim(img_range[0], img_range[1])
    ax.set_ylim(img_range[2], img_range[3])
    ax.axis('off')
    
    if legend:
        ax.add_artist(ax.legend(handles=[obj_dot_1, obj_dot_2, obj_tc_track], 
                                loc=(0.01, 0.01), fontsize=16, framealpha=0.75))
        
    plt.tight_layout()
    return fig, ax

def geo2coord(geo_locs, ref_loc):
    
    # Define UTM projection
    utm_zone = "15N"
    wgs84_proj = "epsg:4326"
    utm_proj = f"epsg:326{utm_zone[:-1]}"
    
    # Convert the reference point to UTM
    tf = Transformer.from_crs(wgs84_proj, utm_proj)
    e0, n0 = tf.transform(ref_loc[1], ref_loc[0])
    
    EN_locs = []
    for i in range(geo_locs.shape[1]):
        loc = geo_locs[:, i]
        e, n = tf.transform(loc[1], loc[0])
        easting = (e - e0) / 1000
        northing = (n - n0) / 1000
        EN_locs.append((easting, northing))

    return np.array(EN_locs).T

def plot_vel_model(df, ax=None, ls='-', lw=2, legend=True):
    
    # Calculate cumulative depth
    depth = df['Thickness (km)'].cumsum()
    depth = np.insert(depth[:-1], [0, len(depth)-1], [0, 2]) * 1e3

    # Properties
    rho = df['rho (g/cm^3)'].to_numpy()
    vp = df['Vp (km/s)'].to_numpy()
    vs = df['Vs (km/s)'].to_numpy()
    rho = np.append(rho, rho[-1])
    vp = np.append(vp, vp[-1])
    vs = np.append(vs, vs[-1])
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(4,4))

    # Plot velocity
    # ax.step(depth, vp, where='post', label='Vp', color='blue', lw=lw, ls=ls)
    ax.step(depth, vs, where='post', label='Vs', color='red', lw=lw, ls=ls)

    # Formatting the plot
    ax.set_xlabel('Depth (m)')
    ax.set_ylabel('Velocity (km/s)')
    ax.grid(True)
    ax.set_xlim([0, 600])
    ax.set_ylim([0, 0.65])

    # Plot effective shear modulus
    mu_bar = rho*vs**2 * (1 - (vs/vp)**2)
    ax1 = ax.twinx()
    ax1.step(depth, mu_bar, where='post', label='$\overline{\mu}$', color='k', lw=lw, ls=ls)
    ax1.set_ylabel('$\overline{\mu}$ (GPa)')
    ax1.set_ylim([0, 0.65])

    # Legend
    if legend:
        handles1, labels1 = ax.get_legend_handles_labels()
        handles2, labels2 = ax1.get_legend_handles_labels()
        handles = handles1 + handles2
        labels = labels1 + labels2
        ax.legend(handles, labels, loc='lower right', fontsize=7)
