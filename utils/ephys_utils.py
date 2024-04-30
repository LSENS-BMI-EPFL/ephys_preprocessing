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
#import neo
#import xarray as xr
#import quantities as pq
#from elephant.conversion import BinnedSpikeTrain

# Modules
sys.path.append(r'C:\Users\bisi\Github\datamanipulation')

def check_if_valid_recording(config, mouse_id, probe_id):
    """
    Check if recording is valid.
    :param config: (dict) config dict.
    :param mouse_id: (str) mouse name.
    :param probe_id: (int) probe id.
    :return:
    """

    path_to_probe_insertion_info = os.path.join(config['mice_info_path'], 'probe_insertion_info.xlsx')
    probe_info_df = pd.read_excel(path_to_probe_insertion_info)
    probe_info = probe_info_df.loc[(probe_info_df['mouse_name'] == mouse_id)
                                      & (probe_info_df['probe_id'] == int(probe_id))]
    if probe_info['valid'].values[0] == 0:
        print('{} probe recording {} not valid. Skipping...'.format(mouse_id, probe_id))
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
    Change stereotaxic coordinates reference for insertion angles.
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

def replace_coil_artefact(spike_array, bin_size, artefact_bins):
    # TODO: to test

    n_neurons = spike_array.shape[0]
    firing_rates = np.nanmean(spike_array / bin_size, axis=1)
    lambdas = firing_rates * bin_size

    poisson_spikes = np.random.poisson(lambdas, size=(n_neurons, len(artefact_bins)))
    spike_array[:, artefact_bins] = poisson_spikes

    return spike_array

def correct_coil_artefact(dataloader, spike_train_array): #TODO: to remove
    """
    Around whisker stimulus time, correct coil artefact by replacing spikes by Poisson spikes.
    :param dataloader: Mouse instance of DataLoader.
    :param spike_train_array: Array of spike trains.
    :return: Corrected spike data array.
    """

    # Define correction window around stimulus time (ms), and baseline window
    stim_time_ms = 1000
    win_pre_ms = 5
    win_post_ms = 6
    correction_start = stim_time_ms - win_pre_ms
    correction_end = stim_time_ms + win_post_ms
    base_fr_start = 0
    base_fr_end = stim_time_ms + win_post_ms

    # Get spike data array, and trial types

    cr_trials = dataloader.get_trial_type_indices(trial_type='CR')
    try:
        whisker_stim_trials = dataloader.datachunk_times[dataloader.datachunk_times['Wh_NoWh'] == 1].index.values
    except KeyError:
        whisker_stim_trials = dataloader.datachunk_times[
            dataloader.datachunk_times['Whisker/NoWhisker'] == 1].index.values

    # Calculate baseline trial-avg. firing rates for each neuron
    spks_for_lambda = np.nanmean(spike_train_array[:, cr_trials, base_fr_start:base_fr_end], axis=1)

    # Generate random Poisson spikes
    rng = np.random.default_rng(seed=None)  # no seed for variability
    poisson_spks = [rng.poisson(lam=spks_for_lambda, size=spks_for_lambda.shape)
                    for i in range(len(whisker_stim_trials))]
    poisson_spks = np.array(poisson_spks).swapaxes(0, 1)

    # Replace spike data in window
    spike_train_array[:, whisker_stim_trials, correction_start:correction_end] = poisson_spks[:, :,
                                                                                 correction_start:correction_end]

    return spike_train_array
