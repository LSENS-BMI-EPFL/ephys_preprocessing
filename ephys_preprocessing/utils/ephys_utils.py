#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: ephys_utils.py
@time: 7/28/2022 10:46 AM
@description: Ephys data manipulation and script utilities.
"""

# Imports
import sys
import re
import os
import pandas as pd
import numpy as np
from pathlib import Path
from collections.abc import Iterable
from loguru import logger
from datetime import datetime

def get_exp_datetime(input_dir):
    """
    Get experiment datetime from input directory.
    :param input_dir:
    :return datetime: str of experiment
    """
    # The session name is the folder before the Ephys folder
    #  e.g. in /scratch/bisi/data/AB133/AB133_20241105_111234/Ephys/, the session name is AB133_20241105_111234

    input_dir = Path(input_dir)

    # Get session folder name
    path_parts = input_dir.parts
    if path_parts[-1].startswith("catgt_"):
        # .../AB133/Recording/AB133_20241105_111234/Ephys/
        session_name = path_parts[-3]
    else:
        # .../AB133/AB133_20241105_111234/Ephys/
        session_name = path_parts[-2]

    # Extract datetime string
    datetime_str = session_name.split('_', 1)[1]

    formatted_date = datetime.strptime(
        datetime_str,
        "%Y%m%d_%H%M%S"
    ).strftime("%d.%m.%Y")
    logger.debug(f"Session {session_name} - Experiment datetime: {formatted_date}")
    return formatted_date



def check_if_valid_recording(config, mouse_id, probe_id, date):
    """
    Check if recording is valid.
    :param config: (dict) config dict.
    :param mouse_id: (str) mouse name.
    :param probe_id: (int) probe id.
    :param date: (str) date of experiment
    :return:
    """
    if mouse_id.startswith('AB') or mouse_id.startswith('MH'):
        path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'joint_probe_insertion_info.xlsx')
        print('Reading experimental metadata:', path_to_probe_insertion_info)
    else:
        path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'probe_insertion_info.xlsx')

    probe_info_df = pd.read_excel(path_to_probe_insertion_info)

    # convert dataframe column
    probe_info_df['date'] = pd.to_datetime(probe_info_df['date'], errors='coerce', dayfirst=True)

    # convert query date
    date_dt = pd.to_datetime(date, dayfirst=True)

    probe_info = probe_info_df.loc[
        (probe_info_df['mouse_name'] == mouse_id) &
        (probe_info_df['probe_id'] == int(probe_id)) &
        (probe_info_df['date'] == date_dt)
     ]

   #probe_info_df['date'] = probe_info_df['date'].astype(str)
   #print(date, type(date), len(probe_info_df[probe_info_df.date==date]))
   #print(probe_info_df.date.values)
   #print(mouse_id,probe_id,date)
   #print(len(probe_info_df.loc[(probe_info_df['mouse_name'] == mouse_id)
   #                                  & (probe_info_df['probe_id'] == int(probe_id))]))
   #probe_info = probe_info_df.loc[(probe_info_df['mouse_name'] == mouse_id)
   #                                  & (probe_info_df['probe_id'] == int(probe_id))
   #                                &  (probe_info_df['date']==date)
   #                                  ]

    # Check if no entries for that mouse
    if probe_info.empty:
        logger.error('No probe insertion info for mouse {} and probe {}. Update probe insertion table.'.format(mouse_id, probe_id))
        return False
    if probe_info['valid'].values[0] == 0:
        logger.warning('Probe insertion for mouse {} and probe {} is not valid. Skipping.'.format(mouse_id, probe_id))
        return False
    return True


def get_probe_version(config,mouse_id,probe_id, date):
    """
    Check if probe_version is 1 or 2.
    :param config:
    :param mouse_id:
    :param probe_id:
    :param date:
    :return:
    """
    if mouse_id.startswith('AB') or mouse_id.startswith('MH'):
        path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'joint_probe_insertion_info.xlsx')
        print('Reading experimental metadata:', path_to_probe_insertion_info)
    else:
        path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'probe_insertion_info.xlsx')

    probe_info_df = pd.read_excel(path_to_probe_insertion_info)

    # convert dataframe column
    probe_info_df['date'] = pd.to_datetime(probe_info_df['date'], errors='coerce', dayfirst=True)

    # convert query date
    date_dt = pd.to_datetime(date, dayfirst=True)

    probe_info = probe_info_df.loc[
        (probe_info_df['mouse_name'] == mouse_id) &
        (probe_info_df['probe_id'] == int(probe_id)) &
        (probe_info_df['date'] == date_dt)
        ]

    probe_info['probe_type'] = probe_info['probe_type'].astype(int)
    # Check if no entries for that mouse
    if probe_info.empty:
        logger.error('No probe insertion info for mouse {} and probe {}. Update probe insertion table.'.format(mouse_id,                                                                                              probe_id))
        return None

    probe_version = probe_info['probe_type'].values[0]
    return probe_version

def convert_stereo_coords(azimuth, elevation):
    """
    Change stereotaxic coordinates reference for insertion angles (for Axel's setup only, AI3209 setup #1)
    Source: https://github.com/petersaj/neuropixels_trajectory_explorer/wiki/General-use
    :param azimuth: (int) azimuth angle as read on L&N setup
    :param elevation: (int) azimuth angle as read on L&N setup
    :return:
    """
    # Convert azimuth angle
    if azimuth < 0:
        azimuth_angle = 360 + azimuth
    elif azimuth > 0:
        azimuth_angle = azimuth
    else: #0°
        azimuth_angle = azimuth

    # Convert elevation angle
    if elevation < 0:
        print('Error, elevation cannot be negative - check probe_insertion sheet.')
        elevation_angle = abs(elevation)
    else:
        elevation_angle = elevation

    return azimuth_angle, elevation_angle

def make_cont_spike_trains(ephys_cluster_df, recording_duration):
    """
    Make list spike trains for a list of clusters for continuous session.
    :param ephys_cluster_df: (pd.DataFrame) from EphysCluster table.
    :param recording_duration: (float) from EphysSession table.
    :return:
    """
    spike_trains_cont = [neo.SpikeTrain(times=c_id.spike_times, t_stop=recording_duration, units='s')
                         for idx, c_id in ephys_cluster_df.iterrows()]
    return spike_trains_cont

def make_binned_cont_spike_trains(spike_trains_cont, bin_size_sec = 0.01):
    """
    Make list of binned spike trains for continuous session.
    :param spike_trains_cont: Output of make_spike_trains.
    :param bin_size_sec: Binning size in seconds.
    :return:
    """
    spike_trains_cont_bin = BinnedSpikeTrain(spike_trains_cont, bin_size= bin_size_sec * pq.s)
    return spike_trains_cont_bin


def make_binned_trial_xarray(spike_trains_cont_bin, trial_outcomes, trial_start_times):
    # Make continuous DataArray
    st_cont_bin_xarr = xr.DataArray(spike_trains_cont_bin.to_array(), dims=('neuron', 'time'))

    trial_xarr_list = []
    bin_size = float(spike_trains_cont_bin.bin_size)
    n_bins = int(1 / bin_size)

    # Slice DataArray at each trial
    for t_start in trial_start_times[:-1]:
        trial_start = int(t_start * n_bins)
        trial_pre = trial_start - n_bins
        trial_post = trial_start + n_bins
        st_trial_bin_xarr = st_cont_bin_xarr.isel(time=slice(trial_pre, trial_post))
        trial_xarr_list.append(st_trial_bin_xarr)

    # Make DataArray
    st_bin_trial_xarr = xr.DataArray(trial_xarr_list, dims=('trial', 'neuron', 'time'),
                                     coords={'trial': trial_outcomes})

    return st_bin_trial_xarr

def flatten_list(l):
    """ Flatten a list of list.
    :param l: A list containing lists.
    :return: Generator of the iterable.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten_list(el)
        else:
            yield el

def extract_ks_version(s):
    """
    Extract Kilosort version number from a string containing 'kilosort' followed by a version number.
    
    Parameters
    ----------
    s : str
        String containing kilosort version (e.g., 'kilosort2', 'kilosort2.5', 'kilosort3')
        
    Returns
    -------
    int or float or None
        Version number as int (e.g., 2, 3) or float (e.g., 2.5) if found, None otherwise
    """
    match = re.search(r'kilosort(\d+(?:\.\d+)?)', s)
    if match:
        version_str = match.group(1)
        return float(version_str) if '.' in version_str else int(version_str)
    return None