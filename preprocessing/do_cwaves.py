#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: do_cwaves.py.py
@time: 30/03/2022 11:40
"""
# Imports
import os
import json
import pprint
import numpy as np
import pandas as pd
import subprocess
import tkinter.filedialog as fdialog
from pathlib import Path
import webbrowser

# Load config
print('Loading C_Waves config:')
with open('cwave_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

# Select raw input ephys recording run
input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory', initialdir=config['data_path'])
session_name = [i for i in os.listdir( os.path.join(input_dir_mouse, 'Recording')) if 'AB' in i][0]
input_dir = os.path.join(input_dir_mouse, 'Recording', session_name, 'Ephys')
print(input_dir)

catgt_epoch_name = os.listdir(input_dir)[0]  # for CatGT processed data
epoch_name = catgt_epoch_name.lstrip('catgt_')

# Select output mouse directory
output_dir_mouse = fdialog.askdirectory(title='Please select mouse output directory', initialdir=config['save_path'])
output_dir = os.path.join(output_dir_mouse, 'Recording', session_name, 'Ephys')
# Get run info and number of probes saved
dirnames = 1

subdir_list = next(os.walk(os.path.join(output_dir, catgt_epoch_name)))[dirnames]
n_probes = len([d for d in subdir_list if 'imec' in d])
print('Run {} recorded using {} probes'.format(epoch_name, n_probes))

# Run C_Waves for each probe
for probe_id in range(n_probes):
    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)

    # Output folder
    path_cwave_output = os.path.join(output_dir, catgt_epoch_name, probe_folder, 'cwaves')
    Path(path_cwave_output).mkdir(parents=True, exist_ok=True)

    # Path to SGLX probe .bin file
    apbin_fname = '{}_tcat.imec{}.ap.bin'.format(epoch_name, probe_id)
    path_to_apbin = os.path.join(input_dir, catgt_epoch_name, probe_folder, apbin_fname)

    # Format cluster files
    path_input_files = os.path.join(output_dir, catgt_epoch_name, probe_folder)

    # Cluster table: incl. peak channel id. Need table row <-> cluster_id equivalence.

    try:
        clus_info = pd.read_csv(os.path.join(path_input_files, 'cluster_info.tsv'), sep='\\t')  # <- assumes Phy performed
    except FileNotFoundError:
        print('Skipping. No spike sorting at', path_input_files)
        continue

    clus_info.set_index(keys='cluster_id', drop=False, inplace=True)  # set index to cluster_id
    clus_table = clus_info.reindex(range(np.max(clus_info.cluster_id) + 1), fill_value=0,
                                   copy=True)  # reindex with missing cluster ids
    clus_table = clus_table[['n_spikes', 'ch']]
    path_clus_table = os.path.join(path_input_files, 'clus_table.npy')
    np.save(path_clus_table, np.array(clus_table.values, dtype=np.int32))

    # Cluster time: spike timestamp (in samples) for each spikes
    spk_times = np.load(os.path.join(path_input_files, 'spike_times.npy'))
    spk_times_df = pd.DataFrame(spk_times, columns=['ts'])  # timestamps col
    path_clus_time = os.path.join(path_input_files, 'clus_time.npy')
    np.save(path_clus_time, spk_times_df['ts'].values)

    # Cluster label: Cluster id for each spike event
    spk_clusters = np.load(os.path.join(path_input_files, 'spike_clusters.npy'))
    spk_clusters_df = pd.DataFrame(spk_clusters, columns=['cluster'])
    path_clus_lbl = os.path.join(path_input_files, 'clus_lbl.npy')
    np.save(path_clus_lbl, spk_clusters_df['cluster'].values)

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
               '-snr_radius={}'.format(config['snr_radius'])
               ]

    print('C_waves command line will run:', command)

    # Run C_waves
    subprocess.run(command, shell=True, cwd=config['base_path'])

    # Open log file
    print('Opening log file')
    webbrowser.open(os.path.join(config['base_path'], 'C_Waves.log'))

    print('Finished C_waves for IMEC probe {} in {}'.format(probe_id, probe_folder))
