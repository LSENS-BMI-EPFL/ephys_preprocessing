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
import os
import pandas as pd
import numpy as np
from collections.abc import Iterable
from loguru import logger


def check_if_valid_recording(config, mouse_id, probe_id, day_id=0):
    """
    Check if recording is valid.
    :param config: (dict) config dict.
    :param mouse_id: (str) mouse name.
    :param probe_id: (int) probe id.
    :param day_index: (int) day index of recordings (naive, expert).
    :return:
    """

    path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'probe_insertion_info.xlsx')
    probe_info_df = pd.read_excel(path_to_probe_insertion_info)
    probe_info = probe_info_df.loc[(probe_info_df['mouse_name'] == mouse_id)
                                      & (probe_info_df['probe_id'] == int(probe_id))
                                      & (probe_info_df['day_of_recording'] == int(day_id))]

    # Check if no entries for that mouse
    if probe_info.empty:
        logger.error('No probe insertion info for mouse {} and probe {}. Update probe insertion table.'.format(mouse_id, probe_id))
        return False
    if probe_info['valid'].values[0] == 0:
        logger.warning('Probe insertion for mouse {} and probe {} is not valid. Skipping.'.format(mouse_id, probe_id))
        return False
    return True

    # Check if there are clusters
    #if cluster_info_df.empty:
    #    return False

    # Check if there are clusters with good or mua labels
    #if not cluster_info_df['group'].isin(['good', 'mua']).any():
    #    return False

    # Check if only noise clusters
    #if cluster_info_df['group'].isin(['noise']).all():
    #    print('Only noise clusters in recording.')
    #    return False

    return


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