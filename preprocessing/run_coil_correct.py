#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_coil_correct.py
@time: 10/31/2023 10:01 AM
"""

# Imports
import os
import pandas as pd
import numpy as np


def main(input_dir):
    """

    :param input_dir:
    :param config:
    :return:
    """

    # Get paths
    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')
    sync_event_folder = os.path.join(input_dir, 'sync_event_times')

    # List IMEC folders
    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Get spike times for recording session
    sync_spike_times_file_list = [f for f in os.listdir(sync_event_folder) if 'spike_times_sec_sync' in f]

    # Get coil times (whisker stimulation) for recording session
    coil_times = np.loadtxt(os.path.join(sync_event_folder, 'whisker_stim_times.txt'))

    # Window of correction in seconds
    pre_stim = 0.0025
    post_stim = 0.0025
    bin_size = pre_stim + post_stim

    # Coil artifact correction for each probe
    for probe_id in probe_ids:
        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        probe_folder_path = os.path.join(input_dir, probe_folder)

        # Get sync spikes times for each probe, in sec
        sync_spikes_times_fname = [f for f in sync_spike_times_file_list if 'imec{}'.format(probe_id) in f][0]
        spike_times = np.load(os.path.join(sync_event_folder, sync_spikes_times_fname))

        # Get spikes times for each cluster
        spk_clusters = np.load(os.path.join(probe_folder_path, 'spike_clusters.npy'))
        spk_clusters_df = pd.DataFrame(spk_clusters, columns=['cluster'])

        # Get firing rates of all clusters
        clus_info_df = pd.read_csv(os.path.join(probe_folder_path, 'cluster_info.tsv'), sep='\\t')
        clus_info_df_corrected = clus_info_df.copy()

        # Iterate over clusters
        for cluster_id in spk_clusters_df.cluster.unique():

            c_fr = clus_info_df['fr'].values[cluster_id]

            # Keep cluster spike times
            c_spike_times = spike_times[spk_clusters_df.cluster == cluster_id]

            # Iterate over stimulation times
            for coil_time in coil_times:

                # Sample spike count from Poisson distribution with mean firing rate
                n_spikes_correct = np.random.poisson(c_fr*bin_size)

                # Shuffle spike times in window, relative to coil time
                spikes_times_sec = np.random.uniform(low=coil_time - pre_stim,
                                                     high=coil_time + post_stim,
                                                     size=n_spikes_correct)

                # Replace spike times with sampled spike times
                c_spike_times_window = c_spike_times[np.where((coil_time - pre_stim < c_spike_times) &
                                                              (c_spike_times < coil_time + post_stim))]

                # Add spikes to window
                window_end_idx = np.where(c_spike_times > coil_time + post_stim)[0]
                c_spike_times = np.insert(c_spike_times, window_end_idx, spikes_times_sec) # this inserts before index

            # Update cluster spike trains and info with all new corrections
            spike_times[spk_clusters_df.cluster == cluster_id] = c_spike_times
            clus_info_df_corrected['n_spikes'].values[cluster_id] = len(c_spike_times)

        # Save corrected spike trains containers
        filename = sync_spikes_times_fname.split('.')[0] + '_coil_corrected.npy'
        np.savetxt(os.path.join(sync_event_folder, filename), spike_times)
        clus_info_df_corrected.to_csv(os.path.join(probe_folder_path, 'cluster_info_coil_corrected.tsv'), sep='\t', index=False)
