#!/usr/bin/env python3

"""
Functions for CM1-related analysis

Author: Qing Ji
"""

# Load python packages
import numpy as np
import netCDF4 as nc
import pandas as pd
from pyproj import Geod
from numpy.polynomial.polynomial import polyfit, polyval


# Create input sounding data for CM1
def create_sounding(df_avp, ignore_wind=True, filename='input_sounding'):
    
    # Extract surface values
    sfc_prs = df_avp['p (mbar)'].iloc[0]
    sfc_theta = df_avp['theta (K)'].iloc[0]
    sfc_wv = df_avp['wv (g/kg)'].iloc[0]

    # Extract the profile data
    pbl_vars = df_avp[['z (m)', 'theta (K)', 'wv (g/kg)', 'u (m/s)', 'v (m/s)']].iloc[1:]
    pbl_vars = pbl_vars.fillna(0)

    # Ignore wind speed
    if ignore_wind:
        pbl_vars[['u (m/s)', 'v (m/s)']] = 0

    # Format the DataFrame
    format_str = f"{{:>10.4f}}"
    pbl_vars = pbl_vars.map(lambda x: format_str.format(x))
    
    # Open the file and write the data
    with open(filename, 'w') as f:
        # Write the surface values
        f.write(f"{sfc_prs:10.4f}\t{sfc_theta:10.4f}\t{sfc_wv:10.4f}\n")
        
        # Write the profile data
        pbl_vars.to_csv(f, sep='\t', header=False, index=False)

    print(f"Sounding data written to {filename}")


# Dry, constant d(theta)/dz sounding (isnd = 8)
# Linearly decreasing wind profile   (iwnd = 8)
def create_ref_sounding(dz=20, p_sfc=1000.0, T_sfc=26.85, lapse=0.005, 
                        rh_sfc=0.90, rh_lapse=-0.001, hurr_V=40, hurr_angle=0):

    # Constants
    g = 9.80665     # m/s^2
    Rd = 287.04     # J/(kg⋅K)
    cp = 1005.7     # J/(kg⋅K)
    p_ref = 1e3     # mbar

    # Scalar height levels (m)
    zh = np.concatenate(([0], np.arange(dz/2, 3000+dz, dz)))

    # Potential temperature (K, lin. decrease with height)
    pi_sfc = (p_sfc / p_ref) ** (Rd/cp)
    th_sfc = (T_sfc + 273.15) / pi_sfc
    theta = th_sfc + lapse * zh

    # Pressure (mbar)
    pi0 = pi_sfc - (g/(cp*lapse)) * np.log(theta/th_sfc)
    pres = p_ref * (pi0**(cp/Rd))

    # Temperature (Celsius)
    T = theta * pi0 - 273.15

    # Relative humidity (N.D., lin. decrease with height)
    rh = rh_sfc + rh_lapse * zh

    # Calculate mixing ratio & specific humidity (g/kg)
    e = e_sw(T) * rh
    wv = 622 * e / (pres*1e2 - e)
    qv = wv / (1 + wv*1e-3)

    # Linearly decreasing wind profile
    wspd = hurr_V * (1.0 - zh/18000.0)
    wspd = wspd * (wspd > 0)

    # Convert to Cartesian components
    hurr_angle_ = np.deg2rad(hurr_angle)
    u = wspd * np.sin(hurr_angle_)
    v = wspd * np.cos(hurr_angle_)

    # Create DataFrame
    df_snd = pd.DataFrame({'z (m)': zh, 'p (mbar)': pres, 
                           'theta (K)': theta, 'wv (g/kg)': wv, 'qv (g/kg)': qv,
                           'u (m/s)': u, 'v (m/s)': v, 'wspd (m/s)': wspd})

    return df_snd


# Saturation vapor pressure (Hardy, 1998)
def e_sw(T):

    # Temperature (K)
    Tk = T + 273.15

    # Series coefficients
    coeff = np.array([-2.8365744e3, -6.028076559e3, 1.954263612e1, 
                      -2.737830188e-2, 1.6261698e-5, 7.0229056e-10, -1.8680009e-13])
    
    # Saturation vapor pressure (Pa)
    es = np.exp(polyval(Tk, coeff) / Tk**2 + 2.7150305 * np.log(Tk))

    return es


# Process QC AVAPS data
def clean_avp(nc_avp):

    # Height levels
    pres = nc_avp.variables['pres'][:]
    inds = np.where(~pres.mask)[0]

    # Geographical location
    lon = nc_avp.variables['lon'][:][inds]
    lat = nc_avp.variables['lat'][:][inds]

    # Read relevant profiles
    pres = pres[inds]                           # mbar
    z = nc_avp.variables['alt'][:][inds]        # m
    T = nc_avp.variables['tdry'][:][inds]       # Celsius
    theta = nc_avp.variables['theta'][:][inds]  # Kelvin
    rh = nc_avp.variables['rh'][:][inds]
    u = nc_avp.variables['u_wind'][:][inds]     # m/s
    v = nc_avp.variables['v_wind'][:][inds]     # m/s
    wspd = nc_avp.variables['wspd'][:][inds]    # m/s

    # Calculate mixing ratio (g/kg)
    e = e_sw(T) * rh / 100
    wv = 622 * e / (pres*1e2 - e)
    qv = wv / (1 + wv*1e-3)

    # Create DataFrame
    df_avp = pd.DataFrame({'z (m)': z, 'lon': lon, 'lat': lat, 'p (mbar)': pres,
                           'theta (K)': theta, 'wv (g/kg)': wv, 'qv (g/kg)': qv,
                           'u (m/s)': u, 'v (m/s)': v, 'wspd (m/s)': wspd})

    return df_avp


# Linear extrapolation of QC AVAPS data
def extrapolate_avp(df_avp, z_top=3000, fit_range=[1500, 2000]):

    # Extract the data to fit
    mask_fit = (df_avp['z (m)'] >= fit_range[0]) & (df_avp['z (m)'] <= fit_range[1])
    z_ = df_avp.loc[mask_fit, 'z (m)']
    theta_ = df_avp.loc[mask_fit, 'theta (K)']
    wv_ = df_avp.loc[mask_fit, 'wv (g/kg)']

    # Linear model
    fit_theta = polyfit(z_, theta_, 1)
    fit_wv = polyfit(z_, wv_, 1)

    # Height levels
    z_cur = df_avp['z (m)'].to_numpy()
    z_new = np.arange(z_cur.max(), z_top, z_cur[-1]-z_cur[-2])  # Add new heights only if necessary
    z = np.concatenate((z_cur, z_new))

    # Extrapolation
    mask_extrap = (z >= fit_range[0])
    z_extrap = z[mask_extrap]
    theta_extrap = polyval(z_extrap, fit_theta).astype('float32')
    wv_extrap = polyval(z_extrap, fit_wv).astype('float32')
    df_extrap = pd.DataFrame({'z (m)': z_extrap, 
                              'theta (K)': theta_extrap, 'wv (g/kg)': wv_extrap})
    
    # Merge data
    df_new = pd.merge(df_avp, df_extrap, on='z (m)', how='outer', suffixes=('', '_extrap'))
    df_new.loc[mask_extrap, 'theta (K)'] = theta_extrap
    df_new.loc[mask_extrap, 'wv (g/kg)'] = wv_extrap
    df_new = df_new.drop(columns=['theta (K)_extrap', 'wv (g/kg)_extrap'])

    return df_new


# Rotate wind speed
def rotate_avp(df_avp, tc_loc):

    # Number of height levels
    Npts = len(df_avp['lon'])

    # Polar coordinate (w.r.t. TC center)
    geod = Geod(ellps='WGS84')
    az, _, _ = geod.inv([tc_loc[0]]*Npts, [tc_loc[1]]*Npts, 
                        df_avp['lon'], df_avp['lat'])
    az_rad = np.deg2rad(az)

    # Radial & Tangential wind
    df_avp_new = df_avp.copy()
    df_avp_new['ur (m/s)'] = df_avp_new['u (m/s)'] * np.sin(az_rad) \
        + df_avp_new['v (m/s)'] * np.cos(az_rad)
    df_avp_new['ut (m/s)'] = -df_avp_new['u (m/s)'] * np.cos(az_rad) \
        + df_avp_new['v (m/s)'] * np.sin(az_rad)
    
    return df_avp_new