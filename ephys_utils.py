#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: ephys_utils.py
@time: 7/28/2022 10:46 AM
"""

# Imports
import sys
import numpy as np

# Modules
sys.path.append(r'C:\Users\bisi\Github\datamanipulation')


def correct_coil_artefact(dataloader, spike_train_array):
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
