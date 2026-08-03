#!/usr/bin/env python3

"""
Functions for reading atmospheric datasets

Author: Qing Ji
"""

# Load python packages
import numpy as np
import pandas as pd

import netCDF4 as nc
from xgrads import CtlDescriptor, open_CtlDataset
import xarray as xr

from . import cm1_analysis as cm1_func


# Read ERA5 reanalysis datasets
def read_era5(era5_file, ts=None, te=None):

    # Select time range
    era5_data = xr.open_dataset(era5_file)
    era5_data = era5_data.sel(valid_time=slice(ts, te))

    # Dimension of dataset
    ndim = len(era5_data.dims)

    # Constants
    R_Earth = 6371229   # m
    g = 9.80665     # m/s^2
    Rd = 287.04     # J/(kg⋅K)
    cp = 1005.7     # J/(kg⋅K)
    p_ref = 1e3     # mbar

    # Single level dataset
    if ndim == 3:
        print('ERA5: Single level')

        # Add wind speed (m/s)
        try:
            s10 = np.sqrt(era5_data['u10']**2 + era5_data['v10']**2)
            era5_data['s10'] = s10
            era5_data['s10'].attrs["units"] = "m/s"
            era5_data['s10'].attrs["long_name"] = "10 metre horizontal wind speed"
        except:
            pass

        # Add relative humidity (N.D.)
        try:
            t2m = era5_data['t2m']
            e2m = cm1_func.e_sw(era5_data['d2m']-273.15)
            RH2m = e2m / cm1_func.e_sw(t2m-273.15)
            era5_data['RH2m'] = RH2m
            era5_data['RH2m'].attrs["units"] = "1"
            era5_data['RH2m'].attrs["long_name"] = "2 metre relative humidity"
        except:
            pass

        # Add elevation (m)
        try:
            gh = era5_data['z'] / g
            alt = R_Earth / (R_Earth - gh) * gh
            era5_data['elev'] = alt
            era5_data['elev'].attrs["units"] = "m"
            era5_data['elev'].attrs["long_name"] = "Elevation"
        except:
            pass

    # Pressure levels dataset
    elif ndim == 4 and ('pressure_level' in era5_data.dims):
        print('ERA5: Pressure level')

        # Add altitude (m)
        gh = era5_data['z'] / g
        alt = R_Earth / (R_Earth - gh) * gh
        era5_data['alt'] = alt
        era5_data['alt'].attrs["units"] = "m"
        era5_data['alt'].attrs["long_name"] = "Altitude"

        # Add wind speed (m/s)
        wspd = np.sqrt(era5_data['u']**2 + era5_data['v']**2)
        era5_data['wspd'] = wspd
        era5_data['wspd'].attrs["units"] = "m/s"
        era5_data['wspd'].attrs["long_name"] = "Horizontal wind speed"

        # Add potential temperature (K)
        theta = era5_data['t'] * (p_ref / era5_data['pressure_level']) ** (Rd / cp)
        era5_data['theta'] = theta
        era5_data['theta'].attrs["units"] = "K"
        era5_data['theta'].attrs["long_name"] = "Potential temperature"

        # Add mixing ratio (g/kg)
        qv = era5_data['q']
        wv = qv / (1 - qv) * 1e3
        era5_data['wv'] = wv
        era5_data['wv'].attrs["units"] = "g/kg"
        era5_data['wv'].attrs["long_name"] = "Mixing ratio"

    return era5_data


# Read MERRA-2 reanalysis datasets
def read_merra(merra_file, ts=None, te=None):

    # Select time range
    merra_data = xr.open_dataset(merra_file)
    merra_data = merra_data.sel(time=slice(ts, te))

    # Dimension of dataset
    ndim = len(merra_data.dims)

    # Constants
    g = 9.80665     # m/s^2
    Rd = 287.04     # J/(kg⋅K)
    cp = 1005.7     # J/(kg⋅K)
    p_ref = 1e3     # mbar

    # Single level dataset
    if ndim == 3:
        print('MERRA-2: Single level')

        # Add wind speed (m/s)
        s10 = np.sqrt(merra_data['U10M']**2 + merra_data['V10M']**2)
        merra_data['s10'] = s10
        merra_data['s10'].attrs["units"] = "m/s"
        merra_data['s10'].attrs["long_name"] = "10 metre horizontal wind speed"

        # Add relative humidity (N.D.)
        q2m = merra_data['QV2M']
        w2m = q2m / (1 - q2m)
        merra_data['wv2m'] = w2m * 1e3
        merra_data['wv2m'].attrs["units"] = "g/kg"
        merra_data['wv2m'].attrs["long_name"] = "2 metre mixing ratio"

    # Pressure levels dataset
    elif ndim == 4 and ('PL' in merra_data.data_vars):
        print('MERRA-2: Model levels')

        # Add wind speed (m/s)
        wspd = np.sqrt(merra_data['U']**2 + merra_data['V']**2)
        merra_data['wspd'] = wspd
        merra_data['wspd'].attrs["units"] = "m/s"
        merra_data['wspd'].attrs["long_name"] = "Horizontal wind speed"

        # Add potential temperature (K)
        theta = merra_data['T'] * (p_ref*1e2 / merra_data['PL']) ** (Rd / cp)
        merra_data['theta'] = theta
        merra_data['theta'].attrs["units"] = "K"
        merra_data['theta'].attrs["long_name"] = "Potential temperature"

        # Add mixing ratio (g/kg)
        qv = merra_data['QV']
        wv = qv / (1 - qv) * 1e3
        merra_data['wv'] = wv
        merra_data['wv'].attrs["units"] = "g/kg"
        merra_data['wv'].attrs["long_name"] = "Mixing ratio"

    elif ndim == 4 and ('PL' not in merra_data.data_vars):
        print('MERRA-2: Pressure levels')

        # Add wind speed (m/s)
        wspd = np.sqrt(merra_data['U']**2 + merra_data['V']**2)
        merra_data['wspd'] = wspd
        merra_data['wspd'].attrs["units"] = "m/s"
        merra_data['wspd'].attrs["long_name"] = "Horizontal wind speed"

        # Add potential temperature (K)
        theta = merra_data['T'] * (p_ref / merra_data['lev']) ** (Rd / cp)
        merra_data['theta'] = theta
        merra_data['theta'].attrs["units"] = "K"
        merra_data['theta'].attrs["long_name"] = "Potential temperature"

        # Add mixing ratio (g/kg)
        qv = merra_data['QV']
        wv = qv / (1 - qv) * 1e3
        merra_data['wv'] = wv
        merra_data['wv'].attrs["units"] = "g/kg"
        merra_data['wv'].attrs["long_name"] = "Mixing ratio"

    return merra_data


# Read radiosonde data
def read_radiosonde(filepath):

    # Constants
    Rd = 287.04     # J/(kg⋅K)
    cp = 1005.7     # J/(kg⋅K)
    p_ref = 1e3     # mbar

    # Read radiosonde CSV file
    df = pd.read_csv(filepath, dtype='str')
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    df = df.replace("", np.nan)
    df = df.astype(float)

    # Drop rows with missing data

    # Extract variables
    z_rs = df['geopotential height_m'].values      # m
    p = df["pressure_hPa"].values                  # hPa
    T = df["temperature_C"].values + 273.15        # K
    r = df["mixing ratio_g/kg"].values             # g/kg
    wspd = df["wind speed_m/s"].values             # m/s

    # Potential temperature (K)
    theta = T * (p_ref/p)**(Rd/cp)

    # Radiosonde profile dictionary
    rs_profile = {'z': z_rs, 'p': p, 'T': T - 273.15, 
                  'wv': r, 'wspd': wspd, 'theta': theta}

    return rs_profile
