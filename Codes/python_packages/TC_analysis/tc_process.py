#!/usr/bin/env python3

"""
Functions for Tropical Cyclone analysis

Author: Qing Ji
"""

# Load python packages
import os, re
import numpy as np
import netCDF4 as nc
import pandas as pd
from pyproj import Geod

from obspy import UTCDateTime
from scipy.interpolate import interp1d

# Tropical cyclone database
base_dir = '/Users/qingji/Documents/Datasets/TC_Track'
tc_dataset = nc.Dataset(os.path.join(base_dir, 'IBTrACS.ALL.v04r00.nc'))


# Search TC in the dataset
def get_tc_info(year, name=None, number=None, agency='USA'):
    if (name is None) and (number is None):
        raise ValueError('Either TC name or number for that year \
                          must be provided!')
    if (name is not None) and (number is not None):
        raise Warning('Receive both TC name and number for that year. \
                       Only searching for TC name!')

    # All TC in the search year
    year_mask = (tc_dataset['season'][:] == year)
    tc_search_list = []
    for obj in tc_dataset['name'][year_mask, :]:
        tc_name = "".join(obj[~obj.mask].astype('str'))
        tc_search_list.append(tc_name)

    # Array index of target TC name
    if (name is not None):
        try:
            tc_ind = np.flatnonzero(year_mask)[tc_search_list.index(name.upper())]
        except Exception:
            raise ValueError('No TC matches the name: %s %d' % (name, year))

    # Array index of target TC number in that year
    else:
        tc_ind = np.flatnonzero(year_mask)[number-1]

    # Create TC info dict
    tc_info = {}

    # ISO time
    time_arr = np.array([])
    for obj in tc_dataset['iso_time'][tc_ind]:
        tc_time = "".join(obj[~obj.mask].astype('str'))
        if tc_time != "":
            time_arr = np.append(time_arr,
                                 UTCDateTime(tc_time, precision=3).datetime)
    tc_info['time'] = time_arr
    max_ind = len(tc_info['time'])

    # TC track
    if agency == 'USA':
        tc_info['lon'] = tc_dataset['usa_lon'][tc_ind][:max_ind]
        tc_info['lat'] = tc_dataset['usa_lat'][tc_ind][:max_ind]

    # Hurricane scale
        tc_info['scale'] = tc_dataset['usa_sshs'][tc_ind][:max_ind]
    # Max sustained wind speed [m/s]
        tc_info['sus_wind'] = tc_dataset['usa_wind'][tc_ind][:max_ind] * 0.514
    # Min sea level pressure [kPa]
        tc_info['pres'] = tc_dataset['usa_pres'][tc_ind][:max_ind] * 0.1

    # Radius of max winds [km] (NOT re-analyzed)
        tc_info['rmw'] = tc_dataset['usa_rmw'][tc_ind][:max_ind] * 1.852
    # Eye diameter [km] (NOT re-analyzed)
        tc_info['eye'] = tc_dataset['usa_eye'][tc_ind][:max_ind] * 1.852

    else:
        tc_info['lon'] = tc_dataset['lon'][tc_ind][:max_ind]
        tc_info['lat'] = tc_dataset['lat'][tc_ind][:max_ind]

        nature_arr = np.array([])
        for obj in tc_dataset['nature'][tc_ind][:max_ind]:
            nature = "".join(obj[~obj.mask].astype('str'))
            nature_arr = np.append(nature_arr, nature)
        tc_info['nature'] = nature_arr

        tc_info['scale'] = tc_dataset['usa_sshs'][tc_ind][:max_ind]
        tc_info['sus_wind'] = tc_dataset['wmo_wind'][tc_ind][:max_ind] * 0.514
        tc_info['pres'] = tc_dataset['wmo_pres'][tc_ind][:max_ind] * 0.1

    # Remove redundant masks
    for key, value in tc_info.items():
        if isinstance(value, np.ma.MaskedArray):
            if not value.mask.any():
                tc_info[key] = value.data

    return tc_info


# Get TC info at specific time
def get_time(tc_info, t):
    if not isinstance(t, UTCDateTime):
        raise TypeError("The input time t must be a UTCDateTime object.")
    t = t.datetime
    
    if t < tc_info['time'].min() or t > tc_info['time'].max():
        raise ValueError("The input time t is outside the TC interval.")
    
    idx = np.where(t >= tc_info['time'])[0][-1]
    return {key: value[[idx, idx+1]] for key, value in tc_info.items()}
    

# Get interpolation function of TC track
def get_interp_func(tc_info):
    
    t0 = tc_info['time'][0]
    t_rel = np.array([(obj-t0).total_seconds() for obj in tc_info['time']])
    f_lon = interp1d(t_rel, tc_info['lon'], bounds_error=False)
    f_lat = interp1d(t_rel, tc_info['lat'], bounds_error=False)
    
    return [f_lon, f_lat], t0


# Interpolate TC track
def interp_track(f_loc, t0, time_range=None, stride=60, time_list=None):
    
    # Convert to UTCDateTime
    t0 = UTCDateTime(t0)

    if (time_range is None) and (time_list is None):
        print(f"Either input time_range or time_list. Both are None now!")
    
    # Timestamp to evaluate TC locations
    if time_list is None:
        tc_reftime = np.arange(0, time_range[1]-time_range[0], 60*stride, dtype=int) \
            + (time_range[0] - t0)
        tc_timestamp = [(t0 + obj).datetime for obj in tc_reftime]
    
    else:
        tc_reftime = np.array([(UTCDateTime(obj) - t0) for obj in time_list])
        tc_timestamp = time_list
    
    # Interpolate TC locations
    tc_lon, tc_lat = f_loc[0](tc_reftime), f_loc[1](tc_reftime)
    
    return tc_timestamp, tc_lon, tc_lat


# Calculate station distance from the hurricane center in km
# Unit in km
def get_station_dist(sta_loc, tc_loc):

    geod = Geod(ellps='WGS84')
    Npts = tc_loc.shape[1]
    _, _, dist = geod.inv([sta_loc[0]]*Npts, [sta_loc[1]]*Npts,
                          tc_loc[0, :], tc_loc[1, :])
    dist = dist / 1e3

    return dist


# Calculate station location under the hurricane center polar coordinate
# Unit in km
def station_geo2polar(sta_info, tc_loc):

    geod = Geod(ellps='WGS84')
    Npts = sta_info.shape[0]
    az, baz, dist = geod.inv([tc_loc[0]]*Npts, [tc_loc[1]]*Npts,
                             sta_info['lon'], sta_info['lat'])
    dist = dist / 1e3

    # Add new columns to station DataFrame
    sta_info_new = sta_info.assign(azimuth=az, dist=dist, backazimuth=baz)
    return sta_info_new


# Read ASOS data
def read_ASOS_data(filename):
    
    asos_raw = pd.read_csv(filename)
    asos_raw.replace(['M','T'], np.nan, inplace=True)
    
    asos_time = np.array(asos_raw['valid'].transform(lambda x: UTCDateTime(x, precision=3).datetime).tolist())
    asos_T = np.array(asos_raw['tmpc'].transform(float).tolist())
    asos_RH = np.array(asos_raw['relh'].transform(float).tolist())
    asos_wind = np.array(asos_raw['sped'].transform(float).tolist()) * 0.44704
    asos_drct = np.array(asos_raw['drct'].transform(float).tolist())
    asos_prs = np.array(asos_raw['alti'].transform(float).tolist()) * 33.8639
    asos_p01 = np.array(asos_raw['p01m'].transform(float).tolist())
    asos_gust = np.array(asos_raw['gust_mph'].transform(float).tolist()) * 0.44704
    asos_theta = (asos_T+273.15) * (1000 / asos_prs) ** 0.286
    
    asos_info = {'time': asos_time, 'T': asos_T, 'theta': asos_theta, 'RH': asos_RH, 'wind': asos_wind, 
                 'drct': asos_drct, 'prs': asos_prs, 'p01': asos_p01, 'gust': asos_gust}
    asos_info = pd.DataFrame(asos_info)

    asos_loc = np.array([asos_raw.iloc[0]['lon'], asos_raw.iloc[0]['lat']])
    
    return asos_info, asos_loc


# Read H*Wind data
def read_Hwind_data(filename):
    
    # Grid data
    df_wind = pd.read_csv(filename, header=3, index_col=0)
    Npts = df_wind.shape[0]

    # Get storm center location
    with open(filename, 'r') as file:
        line = file.readlines()
    line = line[2]
    _lon = re.search(r'(-?\d+\.\d+)', line)
    _lat = re.search(r'(-?\d+\.\d+)', line.split(' and ')[1])

    tc_lon = float(_lon.group()) if _lon else None
    tc_lat = float(_lat.group()) if _lat else None
    
    # Polar coordinate (w.r.t. TC center)
    geod = Geod(ellps='WGS84')
    az, _, dist = geod.inv([tc_lon]*Npts, [tc_lat]*Npts, 
                           df_wind['lon(deg)'], df_wind['lat(deg)'])
    az_rad = np.deg2rad(az)
    
    # Update dataframe
    df_wind['dist(km)'] = dist / 1e3
    df_wind['az(deg)'] = az
    df_wind['R Comp'] = df_wind['U Comp'] * np.sin(az_rad) \
        + df_wind['V Comp'] * np.cos(az_rad)
    df_wind['T Comp'] = -df_wind['U Comp'] * np.cos(az_rad) \
        + df_wind['V Comp'] * np.sin(az_rad)
    
    return df_wind, [tc_lon, tc_lat]


# Obtain H*Wind profile over an azimuth
def get_Hwind_profile(df_wind, az0, az_bin=5, r_bin=1):

    # Azimuthal range
    az_min, az_max = az0 - az_bin/2, az0 + az_bin/2
    az_min = wrap_angle(az_min, method='sym')
    az_max = wrap_angle(az_max, method='sym')
    if az_min < az_max:
        mask = (df_wind['az(deg)'] > az_min) & (df_wind['az(deg)'] < az_max)
    else:
        mask = (df_wind['az(deg)'] > az_min) & (df_wind['az(deg)'] < 180) | \
            (df_wind['az(deg)'] > -180) & (df_wind['az(deg)'] < az_max)
    df_mask = df_wind[mask].copy()

    # Radial bin
    r_bin_edges = np.arange(0, df_mask['dist(km)'].max() + r_bin, r_bin)
    df_mask['r_bin'] = pd.cut(df_mask['dist(km)'], bins=r_bin_edges, right=False)

    # Aggregate within each radial bin
    df_profile = df_mask.groupby('r_bin', observed=False).agg({
        'R Comp': 'mean', 'T Comp': 'mean', 'ws(m/s)': 'mean',
        'dist(km)': 'mean', 'az(deg)': 'mean'
        }).reset_index()

    return df_profile, df_mask


# Wrap angles in degree
def wrap_angle(angle, method='sym'):

    if method == 'sym':
        wrapped_angle = (angle + 180) % 360 - 180
        return wrapped_angle
    
    if method == 'asym':
        wrapped_angle = angle % 360
        return wrapped_angle
