#!/usr/bin/env python3

"""
Functions for downloading seismic data

Author: Qing Ji
"""

# Load python packages
import os, logging
import numpy as np
import fnmatch

from obspy.core.inventory import PolynomialResponseStage
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.mass_downloader import Restrictions, MassDownloader, GlobalDomain


# Download seismic, pressure, etc. record (one trace)
# Instrument response can be removed if ordered
# Pressure record in Pa or mbar, Displacement record in μm (positive for downward)
def download_trace(sta_char, time_range, channel='LHZ', resp_to='DISP',
                   pre_filt='default', extra_portion=0.1, client=Client("IRIS"),
                   prs_unit='Pa', rm_resp=True):

    if pre_filt == 'default':
        pre_filt = [5e-4, 1e-3, 45, 50]

    # Time range
    ts = time_range[0]
    te = time_range[1]

    # Station codes
    net = sta_char[0]
    stnm = sta_char[1]
    loc = '*' if len(sta_char) < 3 else sta_char[2]

    # Obtain station location and polarity
    sta_inv = client.get_stations(starttime=ts, endtime=te, level='channel',
                                  network=sta_char[0], station=sta_char[1],
                                  channel=channel)
    cha_info = sta_inv[0][0][0]
    sta_loc = np.array([cha_info.longitude, cha_info.latitude])

    ### Not removing instrument response, just return raw data
    if not rm_resp:
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.stats.unit = 'Counts'
        return trace, sta_loc

    ### Remove instrument response
    # Infrasound pressure data
    if fnmatch.fnmatch(channel, '*D[FH]'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_response(pre_filt=pre_filt)
        trace.stats.unit = 'Pressure (Pa)'

    # Barometric pressure data
    elif fnmatch.fnmatch(channel, '*DO'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        
        try:
            trace.remove_response(pre_filt=None)
        except:
            # LDO channel with location code '30' or '31'
                # Stage 1: PolynomialResponseStage from Pa to counts, gain: UNKNOWN
                # Stage 2: ResponseStage from counts to counts, gain: 1
                # Stage 3: CoefficientsTypeResponseStage from counts to counts, gain: 1

            # Manually remove PolynomialResponseStage
            response = trace.stats.response
            if isinstance(response.response_stages[0], PolynomialResponseStage):
                if response.response_stages[0].stage_gain is None:
                    gain = 1
                else:
                    gain = response.response_stages[0].stage_gain
                coefficients = response.response_stages[0].coefficients[:]
                for i in range(len(coefficients)):
                    coefficients[i] /= np.power(gain, i)
                trace.data = np.poly1d(coefficients[::-1])(trace.data)

        if prs_unit == 'Pa':
            trace.stats.unit = 'Pressure (Pa)'
        else:
            trace.data = trace.data / 1e2
            trace.stats.unit = 'Pressure (mbar)'

    elif fnmatch.fnmatch(channel, '*DM'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        
        if prs_unit == 'Pa':
            trace.stats.unit = 'Pressure (Pa)'
        else:
            trace.data = trace.data / 1e2
            trace.stats.unit = 'Pressure (mbar)'

    # Seismic data
    elif fnmatch.fnmatch(channel, '*H[EN12Z]'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_response(output='VEL', pre_filt=pre_filt)

        if resp_to == 'DISP':
            # Integration from velocity to displacement
            trace.integrate()
            trace.stats.unit = 'Displacement ($\mu$m)'
        elif resp_to == 'ACC':
            # Differentiate from velocity to acceleration
            trace.differentiate()
            trace.stats.unit = 'Acceleration ($\mu$m/s$^2$)'
        elif resp_to == 'VEL':
            trace.stats.unit = 'Acceleration ($\mu$m/s)'
        trace.detrend()
        
        # Conversion to μm and flip polarity (positive downward)
        trace.stats.dip = cha_info.dip
        if fnmatch.fnmatch(channel, '*HZ'):
            if trace.stats.dip == -90.0:
                trace.data = -trace.data * 1e6
            elif trace.stats.dip == 90.0:
                trace.data = trace.data * 1e6
            else:
                raise ValueError('Weird dip of the vertical channel: %.1f'
                                 % trace.stats.dip)
        
        if fnmatch.fnmatch(channel, '*H[EN12]'):
            trace.data = trace.data * 1e6

    # Wind speed data
    elif fnmatch.fnmatch(channel, '*WS'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Wind Speed (m/s)'

    # Wind direction data
    elif fnmatch.fnmatch(channel, '*WD'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Wind Direction (deg)'

    # Temperature data
    elif fnmatch.fnmatch(channel, '*IO'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Humidity (%)'

    # Humidity data
    elif fnmatch.fnmatch(channel, '*KO'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Temperature (°C)'

    # Rainfall data
    elif fnmatch.fnmatch(channel, '*RO'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Rainfall (mm/hr)'

    # Hail data
    elif fnmatch.fnmatch(channel, '*RH'):
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        trace.remove_sensitivity()
        trace.stats.unit = 'Hail (#/cm$^2$/hr)'


    else:
        trace = client.get_waveforms(net, stnm, loc, channel,
                                     attach_response=True,
                                     starttime=ts-(te-ts)*extra_portion,
                                     endtime=te+(te-ts)*extra_portion)
        trace = trace.merge()[0]
        # trace.remove_response(output='VEL', pre_filt=pre_filt)

        # Integration from velocity to displacement
        # trace.integrate()
        # trace.detrend()
        trace.stats.unit = 'Displacement ($\mu$m)'
        trace.data = trace.data * 1e-6

    return trace, sta_loc


# Check number of data points
def check_npts(stream):
    
    npts_list = [trace.stats.npts for trace in stream]
    
    if len(set(npts_list)) > 1:
        min_npts = min(npts_list)
        
        for trace, npts in zip(stream, npts_list):
            if npts == min_npts:
                pass
            
            # Infrasound LDF channel sometimes have 1 data point less than
            # seismic channels LH*
            elif npts == min_npts + 1:
                trace.data = trace.data[1:]
                trace.stats.starttime += trace.stats.delta
            
            else:
                raise Exception('npts differ more than 1 among traces!')

    return stream


### Mass downloading functions ###

# Directory to store waveforms
data_dir = '/Users/qingji/Documents/Datasets/tmp/waveform'
resp_dir = '/Users/qingji/Documents/Datasets/tmp/response'

# Mass downloader
def mdl_data(sta_char, ts, te, channel, minimum_length=0.5, client=None,
             chunk_length=86400):
    """
    Download data for one specific station

    Parameters:
        sta_char (list): List with format [network, station]
        ts, te (str): Time range to download data
        channel_list (list): List of channels to download data
        station_time (list): List of UTCDateTime for the start and end time 
                             of the station xml file
        chunk_length (int or None): Length of each chunk in seconds. None: No chunk

    No returns. Check data under ./{tmp_dir}/{network}.{station}/
    """
    
    # Mass downloader
    domain = GlobalDomain()
    if client is None:
        mdl = MassDownloader(providers=["IRIS"])
    else:
        mdl = MassDownloader(providers=client)

    # Logger
    logger = logging.getLogger("obspy.clients.fdsn.mass_downloader")
    logger.setLevel(logging.ERROR)
    for handler in logger.handlers:
        handler.setLevel(logging.ERROR)
    
    # Download MSEED data
    if chunk_length is None:
        restrictions = Restrictions(starttime=ts, endtime=te,
                                    network=sta_char[0], station=sta_char[1], 
                                    channel=channel, chunklength_in_sec=None,
                                    reject_channels_with_gaps=False,
                                    sanitize=False, minimum_length=0.0,
                                    location_priorities=['EP','','00','10'],
                                    minimum_interstation_distance_in_m=0)
        mdl.download(domain, restrictions, mseed_storage=get_mseed_storage_full,
                     stationxml_storage=get_xml_path)

    else:
        restrictions = Restrictions(starttime=ts, endtime=te,
                                    network=sta_char[0], station=sta_char[1], 
                                    channel=channel, chunklength_in_sec=86400,
                                    reject_channels_with_gaps=False,
                                    sanitize=False, minimum_length=minimum_length,
                                    location_priorities=['EP','','00','10'],
                                    minimum_interstation_distance_in_m=0)
        mdl.download(domain, restrictions, mseed_storage=get_mseed_storage_julian,
                     stationxml_storage=get_xml_path)


# Function: Get directory path for waveform data
def get_mseed_dir(network, station, 
                  location=None, channel=None, starttime=None, endtime=None):
    """
    Custom function for ObsPy mass downloader.

    Returns:
        mseed_dir (str): Directory path to save waveform data
    """
    
    # File path for saving data
    mseed_dir = os.path.join(data_dir, ".".join([network, station]))
    
    return mseed_dir


# Function: Get file path for waveform data
def get_mseed_storage(network, station, location, channel, starttime, endtime):
    """
    Custom function for ObsPy mass downloader.

    Returns:
        mseed_file (str/bool): File path to be saved as. If True, then the MiniSEED 
        file is assumed to already be available and will not be downloaded again
    """
    
    # File path for saving data
    ts_mod = starttime + 0.5
    mseed_file = os.path.join(get_mseed_dir(network, station), channel,
                              "%4d%02d%02d.%s.mseed" %(ts_mod.year, ts_mod.month, ts_mod.day, channel))
    
    # If file already exists, do not download again
    if os.path.exists(mseed_file):
        return True
    
    return mseed_file


# Function: Get file path for waveform data (julian day filename)
def get_mseed_storage_julian(network, station, location, channel, starttime, endtime):
    """
    Custom function for ObsPy mass downloader.

    Returns:
        mseed_file (str/bool): File path to be saved as. If True, then the MiniSEED 
        file is assumed to already be available and will not be downloaded again
    """
    
    # File path for saving data
    ts_mod = starttime + 0.5
    mseed_file = os.path.join(get_mseed_dir(network, station), channel,
                              "%4d%03d.%s.mseed" %(ts_mod.year, ts_mod.julday, channel))
    
    # If file already exists, do not download again
    if os.path.exists(mseed_file):
        return True
    
    return mseed_file


# Function: Get file path for waveform data (julian day filename)
def get_mseed_storage_full(network, station, location, channel, starttime, endtime):
    """
    Custom function for ObsPy mass downloader.

    Returns:
        mseed_file (str/bool): File path to be saved as. If True, then the MiniSEED 
        file is assumed to already be available and will not be downloaded again
    """
    
    # File path for saving data
    mseed_file = os.path.join(get_mseed_dir(network, station),
                              "%s.%s.%s.mseed" %(network, station, channel))
    
    # If file already exists, do not download again
    if os.path.exists(mseed_file):
        return True
    
    return mseed_file


# Function: Get file path for station xml data
def get_xml_path(network, station, channels, starttime, endtime):
    """
    Custom function for ObsPy mass downloader.

    Returns:
        xml_path (str): File path to save station xml file
    """
    
    # Station channels
    available_channels = []
    missing_channels = []
    channel_str = ''
    for location, channel in channels:
        missing_channels.append((location, channel))
        channel_str = channel_str + channel

    # File path for saving station xml data
    xml_path = os.path.join(resp_dir, channel_str, f"{network}.{station}.xml") 
    
    return {"available_channels": available_channels, 
            "missing_channels": missing_channels, 
            "filename": xml_path}