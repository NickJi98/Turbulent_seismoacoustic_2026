#!/usr/bin/env python3

"""
Plotting functions for Tropical Cyclone analysis

Author: Qing Ji
"""

# Load python packages
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from pyproj import Geod

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

from seismic.seis_station import get_sta_info

# Obspy client
client_default = Client("IRIS")


# Plot TC track (and seismic stations)
def plot_tc_track(tc_info, extent=None, category=False, sta_info=None,
                  client=None, plot_sta=False, stnm_list=None, aspect=None):

    # IRIS client
    if client is None:
        client = client_default

    # Map region
    if extent is None:
        margin_in_deg = [10, 5]
        extent = [np.min(tc_info['lon'])-margin_in_deg[0],
                  np.max(tc_info['lon'])+margin_in_deg[0],
                  np.min(tc_info['lat'])-margin_in_deg[1],
                  np.max(tc_info['lat'])+margin_in_deg[1]]

    # Map projection
    proj = ccrs.Mercator()
    proj_data = ccrs.PlateCarree()

    fig = plt.figure()
    ax = plt.axes(projection=proj)
    ax.set_extent(extent, crs=proj_data)
    ax.add_feature(cfeature.STATES, edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    # Plot TC track
    ax.plot(tc_info['lon'], tc_info['lat'], 'b-', transform=proj_data)

    if category is False:
        ax.scatter(tc_info['lon'], tc_info['lat'], s=60, c='blue', marker='o',
                   transform=proj_data, label='Hurricane Track')
    else:
        ax.scatter(tc_info['lon'][tc_info['scale'] >= 1],
                   tc_info['lat'][tc_info['scale'] >= 1],
                   s=60, c='orange', marker='o',
                   transform=proj_data, label='Hurricane')
        ax.scatter(tc_info['lon'][tc_info['scale'] == 0],
                   tc_info['lat'][tc_info['scale'] == 0],
                   s=60, c='green', marker='o',
                   transform=proj_data, label='Tropical Storm')
        ax.scatter(tc_info['lon'][tc_info['scale'] == -1],
                   tc_info['lat'][tc_info['scale'] == -1],
                   s=60, c='blue', marker='o',
                   transform=proj_data, label='Tropical Depression')
        ax.scatter(tc_info['lon'][tc_info['scale'] < -1],
                   tc_info['lat'][tc_info['scale'] < -1],
                   s=60, c='gray', marker='o', transform=proj_data)

    # Plot stations
    if plot_sta:
        if sta_info is None:
            print('Download seismic station information ...')
            sta_inv = client.get_stations(starttime=UTCDateTime(tc_info['time'][0]),
                                          endtime=UTCDateTime(tc_info['time'][-1]),
                                          minlongitude=extent[0],
                                          maxlongitude=extent[1],
                                          minlatitude=extent[2],
                                          maxlatitude=extent[3],
                                          channel='LHZ,BHZ',
                                          level='channel')
            sta_info = get_sta_info(sta_inv)
            print('Finish.')

        ax.scatter(sta_info['lon'], sta_info['lat'],
                   s=60, c='k', marker='^', edgecolors='none',
                   transform=proj_data, label='Station')

        # Plot all station names
        if stnm_list == 'all':
            transform = proj_data._as_mpl_transform(ax)
            for i, stnm in enumerate(sta_info['stnm']):
                if (sta_info['lon'][i] < extent[0]) \
                        or (sta_info['lon'][i] > extent[1]) \
                        or (sta_info['lat'][i] < extent[2]) \
                        or (sta_info['lat'][i] > extent[3]):
                    continue
                ax.annotate(stnm, (sta_info['lon'][i], sta_info['lat'][i]),
                            xycoords=transform)
        
        # Plot specific station names
        elif stnm_list is not None:
            for stnm_ in stnm_list:
                one_info = sta_info[sta_info['stnm']==stnm_]
                if one_info.empty:
                    print(f"Station {stnm_} is not found.")
                    continue
                ax.annotate(stnm_, (one_info.iloc[0]['lon'], one_info.iloc[0]['lat']),
                            xycoords=proj_data, fontsize=14)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--')
    gl.xlabel_style = {'size': 18}
    gl.ylabel_style = {'size': 18}
    gl.top_labels = False
    gl.right_labels = False
    ax.legend()

    if aspect is None:
        aspect = 1.0

    ax.set_aspect(aspect)
    fig.set_size_inches(10/aspect, 10)

    return [fig, ax], sta_info


# Plot station distribution
def plot_stations(sta_info, extent=None, plot_stnm=False, aspect=None,
                  markersize=60, fontsize=18):

    # Map region
    if extent is None:
        margin_in_deg = [10, 5]
        extent = [np.min(sta_info['lon'])-margin_in_deg[0],
                  np.max(sta_info['lon'])+margin_in_deg[0],
                  np.min(sta_info['lat'])-margin_in_deg[1],
                  np.max(sta_info['lat'])+margin_in_deg[1]]

    # Map projection
    proj = ccrs.Mercator()
    proj_data = ccrs.PlateCarree()

    fig = plt.figure()
    ax = plt.axes(projection=proj)
    ax.set_extent(extent, crs=proj_data)
    ax.add_feature(cfeature.STATES, edgecolor='gray')
    ax.add_feature(cfeature.COASTLINE)

    ax.scatter(sta_info['lon'], sta_info['lat'],
               s=markersize, c='k', marker='^', edgecolors='none',
               transform=proj_data, label='Station')

    if plot_stnm:
        transform = proj_data._as_mpl_transform(ax)
        sta_info = sta_info.reset_index()
        for i, stnm in enumerate(sta_info['stnm']):
            if (sta_info['lon'][i] < extent[0]) \
                    or (sta_info['lon'][i] > extent[1]) \
                    or (sta_info['lat'][i] < extent[2]) \
                    or (sta_info['lat'][i] > extent[3]):
                continue
            ax.annotate(stnm, (sta_info['lon'][i], sta_info['lat'][i]),
                        xycoords=transform, fontsize=fontsize)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--')
    gl.xlabel_style = {'size': 18}
    gl.ylabel_style = {'size': 18}
    gl.top_labels = False
    gl.right_labels = False
    ax.legend()

    if aspect is None:
        aspect = 1.0

    ax.set_aspect(aspect)
    fig.set_size_inches(10/aspect, 10)

    return [fig, ax]


# Plot TC translation speed
def plot_tc_trans(tc_info):

    # TC track segment length
    geod = Geod(ellps='WGS84')
    _, _, seg_dist = geod.inv(tc_info['lon'][:-1], tc_info['lat'][:-1],
                              tc_info['lon'][1:], tc_info['lat'][1:])

    # Mid-time and time interval
    dt = tc_info['time'][1:] - tc_info['time'][:-1]
    dt_int = [obj.total_seconds() for obj in dt]
    timestamp = tc_info['time'][1:] + dt/2

    # Translation speed in m/s
    speed_trans = seg_dist / dt_int

    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(timestamp, speed_trans, 'k-')
    ax.xaxis_date()
    ax.set_ylabel('Speed [m/s]')
    ax.set_ylim(0, np.max(speed_trans))
    ax.set_title('Hurricane Translation Speed')
    ax.grid()
    fig.show()

    return timestamp, speed_trans


# Plot station distance from the hurricane center
def plot_station_dist(timestamp, dist, title='Title'):

    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(timestamp, dist, 'k-')
    ax.xaxis_date()
    ax.set_ylabel('Distance [km]')
    ax.set_ylim(0, np.nanmax(dist))
    ax.set_title('Distance from Hurricane Center, %s' % title)
    ax.grid()
    fig.show()


# Plot HWind snapshot
def plot_Hwind_snapshot(df_wind):
    
    # Reshape data
    Nx = np.sqrt(df_wind.shape[0]).astype(int)
    lon_mat = df_wind['lon(deg)'].to_numpy().reshape(Nx, Nx)
    lat_mat = df_wind['lat(deg)'].to_numpy().reshape(Nx, Nx)
    ws = df_wind['ws(m/s)'].to_numpy().reshape(Nx, Nx)
    
    fig, ax = plt.subplots(figsize=(4,3))
    obj = ax.pcolormesh(lon_mat, lat_mat, ws, cmap='jet', shading='gouraud')
    
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    
    fig.colorbar(obj, ax=ax, label='Wind Speed (m/s)')
    
    return fig, ax