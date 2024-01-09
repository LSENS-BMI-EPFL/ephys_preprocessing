#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_cwaves.py
@time: 8/2/2023 11:13 PM
"""

# Imports
import os
import pathlib
import subprocess
import webbrowser
import pandas as pd
import numpy as np


def main(input_dir, config):
    """

    :param input_dir:
    :param config:
    :return:
    """

    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Run C_Waves for each probe
    for probe_id in probe_ids:

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)

        # Create output folder
        path_cwave_output = os.path.join(input_dir, probe_folder, 'cwaves')
        pathlib.Path(path_cwave_output).mkdir(parents=True, exist_ok=True)

        # Path to  probe binary file
        apbin_fname = '{}_tcat.imec{}.ap.bin'.format(epoch_name, probe_id)
        path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)

        # Prepare spike data for C_Waves
        path_input_files = os.path.join(input_dir, probe_folder)

        # Create cluster table: incl. peak channel id. Need table row number <-> cluster_id equivalence (hence reindexing below).
        try:
            clus_info = pd.read_csv(os.path.join(path_input_files, 'cluster_info.tsv'),   # <- assumes Phy performed
                                    sep='\\t')
        except FileNotFoundError:
            print('Skipping probe. No spike sorting at', path_input_files)
            continue

        clus_info.set_index(keys='cluster_id', drop=False, inplace=True)  # set index to cluster_id
        clus_table = clus_info.reindex(range(np.max(clus_info.cluster_id) + 1),
                                       fill_value=0,
                                       copy=True)  # reindex with missing cluster ids
        clus_table = clus_table[['n_spikes', 'ch']]
        path_clus_table = os.path.join(path_input_files, 'clus_table.npy')
        np.save(path_clus_table, np.array(clus_table.values, dtype=np.uint32))

        # Cluster time: spike timestamp (in samples) for each spikes
        spk_times = np.load(os.path.join(path_input_files, 'spike_times.npy'))
        spk_times_df = pd.DataFrame(spk_times, columns=['ts'])  # timestamps col
        path_clus_time = os.path.join(path_input_files, 'clus_time.npy')
        np.save(path_clus_time, spk_times_df['ts'].values)

        # Cluster label: Cluster id for each spike event
        spk_clusters = np.load(os.path.join(path_input_files, 'spike_clusters.npy'))
        spk_clusters_df = pd.DataFrame(spk_clusters, columns=['cluster'])
        path_clus_lbl = os.path.join(path_input_files, 'clus_lbl.npy')
        np.save(path_clus_lbl, np.array(spk_clusters_df['cluster'].values, dtype=np.uint32))


        # Write C_Waves command
        command = ['C_waves',
                   '-spikeglx_bin={}'.format(path_to_apbin),
                   '-clus_table_npy={}'.format(path_clus_table),
                   '-clus_time_npy={}'.format(path_clus_time),
                   '-clus_lbl_npy={}'.format(path_clus_lbl),
                   '-dest={}'.format(path_cwave_output),
                   '-samples_per_spike={}'.format(config['samples_per_spike']),
                   '-pre_samples={}'.format(config['pre_samples']),
                   '-num_spikes={}'.format(config['num_spikes']),
                   '-snr_radius_um={}'.format(config['snr_radius']) # requires snsGeomMap metadata entry, otherwise use -sns_radius
                   ]

        print('C_waves command line will run:', command)

        print('Running C_waves for IMEC probe {}...'.format(probe_id))
        subprocess.run(command, shell=True, cwd=config['cwaves_path'])

        print('Opening log file')
        webbrowser.open(os.path.join(config['cwaves_path'], 'C_Waves.log'))
        # Remove useless mean_waveform rows (necessary for C_Waves to run), then resave
        mean_waveforms = np.load(os.path.join(path_cwave_output, 'mean_waveforms.npy'))
        clus_table['temp_matching'] = clus_table.apply(lambda x: x.n_spikes == x.ch == 0, axis=1)
        ids_to_remove = clus_table[clus_table['temp_matching'] == True].index
        mean_waveforms = np.delete(mean_waveforms, ids_to_remove, axis=0)

        np.save(os.path.join(path_cwave_output, 'mean_waveforms.npy'), mean_waveforms)

    return
