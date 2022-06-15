#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
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
from collections.abc import Iterable
import webbrowser


# Load config
print('Loading C_Waves config')
with open('cwave_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

# Select raw input ephys recording run
input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory', initialdir=config['data_path'])
input_dir = os.path.join(input_dir_mouse, 'Recording/Ephys')

epoch_name = os.listdir(input_dir)[0] #for raw data
catgt_epoch_name = 'catgt_{}'.format(epoch_name) #for CatGT processed data
run_name = os.listdir(input_dir)[0][0:-3] #mouse name

# Select output mouse directory
output_dir_mouse = fdialog.askdirectory(title='Please select mouse output directory', initialdir=config['save_path'])
output_dir = os.path.join(output_dir_mouse, 'Recording/Ephys')

# Get run info and number of probes saved
dirnames = 1
n_probes = len(next(os.walk(os.path.join(output_dir, catgt_epoch_name)))[dirnames])
print('Run {} recorded using {} probes'.format(epoch_name, n_probes))

# Run C_Waves for each probe
for probe_id in range(n_probes):
    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)

    #Output folder
    path_cwave_output = os.path.join(output_dir, catgt_epoch_name, probe_folder, 'cwaves')
    Path(path_cwave_output).mkdir(parents=True, exist_ok=True)

    #Path to SGLX probe .bin file
    apbin_fname = '{}_t0.imec{}.ap.bin'.format(epoch_name, probe_id)
    path_to_apbin = os.path.join(input_dir, epoch_name, probe_folder, apbin_fname)

    #Format cluster files
    path_input_files = os.path.join(output_dir, catgt_epoch_name, probe_folder, 'ks25')

    #Cluster table: incl. peak channel id. Need table row <-> cluster_id equivalence.
    clus_info = pd.read_csv(os.path.join(path_input_files,'cluster_info.tsv'), sep='\\t') #<- assumes Phy performed
    clus_info.set_index(keys='cluster_id', drop=False, inplace=True) #set index to cluster_id
    clus_table = clus_info.reindex(range(np.max(clus_info.cluster_id) + 1), fill_value=0, copy=True) #reindex with missing cluster ids
    clus_table = clus_table[['n_spikes', 'ch']]
    path_clus_table = os.path.join(path_input_files, 'clus_table.npy')
    np.save(path_clus_table, np.array(clus_table.values, dtype=np.int32))

    #Cluster time: spike timestamp (in samples) for each spikes
    spk_times = np.load(os.path.join(path_input_files, 'spike_times.npy'))
    spk_times_df = pd.DataFrame(spk_times, columns=['ts'])  # timestamps col
    path_clus_time = os.path.join(path_input_files, 'clus_time.npy')
    np.save(path_clus_time, spk_times_df['ts'].values)

    #Cluster label: Cluster id for each spike event
    spk_clusters = np.load(os.path.join(path_input_files, 'spike_clusters.npy'))
    spk_clusters_df = pd.DataFrame(spk_clusters, columns=['cluster'])
    path_clus_lbl = os.path.join(path_input_files, 'clus_lbl.npy')
    np.save(path_clus_lbl, spk_clusters_df['cluster'].values)

    #Write C_Waves command
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
