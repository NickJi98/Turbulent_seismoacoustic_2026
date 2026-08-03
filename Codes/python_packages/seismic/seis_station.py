#!/usr/bin/env python3

"""
Functions for organizing station metadata

Author: Qing Ji
"""

# Load python packages
import numpy as np
import pandas as pd

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

# Obspy client
client_default = Client("IRIS")


# Obtain station location DataFrame based on Inventory
def get_sta_loc(sta_inv):

    # Lists that store information
    net_list, stnm_list = [], []
    lat_list, lon_list = [], []
    elev_list, cha_list = [], []

    for network in sta_inv:
        for station in network:
            net_list.append(network.code)
            stnm_list.append(station.code)
            lat_list.append(station.latitude)
            lon_list.append(station.longitude)
            elev_list.append(station.elevation)
            cha_list.append('.'.join([cha.code for cha in station]))

    # Create dataframe
    sta_loc = pd.DataFrame({'stnm': stnm_list, 'net': net_list,
                            'lon': np.array(lon_list), 'lat': np.array(lat_list),
                            'elev': np.array(elev_list), 'cha': cha_list})
    return sta_loc


# Obtain station DataFrame based on Inventory
def get_sta_info(sta_inv):

    # Lists that store information
    net_list, stnm_list = [], []
    lat_list, lon_list = [], []
    elev_list, dep_list = [], []
    cha_list = []

    for network in sta_inv:
        net = network.code
        for station in network:
            sta = station.code
            for channel in station:
                net_list.append(net)
                stnm_list.append(sta)
                lat_list.append(channel.latitude)
                lon_list.append(channel.longitude)
                elev_list.append(channel.elevation)
                dep_list.append(channel.depth)
                cha_list.append(channel.code)

    # Create dataframe
    sta_info = pd.DataFrame({'stnm': stnm_list, 'net': net_list,
                             'lon': np.array(lon_list), 'lat': np.array(lat_list),
                             'elev': np.array(elev_list), 'depth': np.array(dep_list),
                             'channel':np.array(cha_list)})
    return sta_info


# Pick out stations having co-located channels
def merge_sta_info(sta_info1, sta_info2):

    sta_list = sta_info2['stnm'].tolist()
    ind_list = []

    for ind in range(len(sta_info1)):
        sta = sta_info1.loc[ind, 'stnm']
        if sta in sta_list:
            ind_list.append(ind)

    sta_info_merge = sta_info1.iloc[ind_list, :4]
    sta_info_merge = sta_info_merge.reset_index(drop=True)

    return sta_info_merge


# Search available stations
def search_station(loc, domain='rect', dist=[10,10], ts=None, te=None, 
                   network=None, station=None, location=None, channel='LHZ,LDF', 
                   level='response', client=None):

    # IRIS client
    if client is None:
        client = client_default

    # Center location
    lon_c, lat_c = loc[0], loc[1]

    # Search domain (distance in degrees)
    if domain == 'rect':
        minlon, maxlon = lon_c-dist[0]/2, lon_c+dist[0]/2
        minlat, maxlat = max(lat_c-dist[1]/2, -90), min(lat_c+dist[1]/2, 90)
        minlon = (minlon + 180) % 360 - 180
        maxlon = (maxlon + 180) % 360 - 180

        sta_inv = client.get_stations(starttime=ts, endtime=te, network=network,
                                      station=station, location=location, channel=channel,
                                      level=level, minlongitude=minlon, maxlongitude=maxlon,
                                      minlatitude=minlat, maxlatitude=maxlat)

    elif domain == 'circ':
        sta_inv = client.get_stations(starttime=ts, endtime=te, network=network,
                                      station=station, location=location, channel=channel,
                                      level=level, longitude=loc[0], latitude=loc[1], 
                                      minradius=dist[0], maxradius=dist[1])
        
    else:
        raise ValueError("Argument 'domain' should be either 'rect' or 'circ'!")
    
    return sta_inv



