#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_artifact_correction.py
@time: 2/14/2024 11:32 AM
"""


# # Imports
import os
import shutil
import subprocess
import webbrowser

import numpy as np
import pathlib

from matplotlib import pyplot as plt

# Import readers
import readSGLX


def main(input_dir, config):
    """
    Run artifact correction on CatGT-processed ephys data using TPrime-aligned artifact times.
    :param input_dir:
    :param config:
    :return:
    """

    epoch_name = os.listdir(input_dir)[0]
    run_name = epoch_name[6:]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    n_probes = len(probe_folders)

    print('Running artifact correction...')
    for probe_id in range(n_probes):
        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)

        # Get ap-band binary data
        ap_bin_filename = [f for f in os.listdir(probe_path) if 'ap.bin' in f][0]
        ap_meta_filename = [f for f in os.listdir(probe_path) if 'ap.meta' in f][0]
        ap_bin_path = pathlib.Path(probe_path, ap_bin_filename)
        ap_meta_dict = readSGLX.readMeta(pathlib.Path(probe_path, ap_meta_filename))

        # Create a copy of the meta file with the corrected name (needed by SpikeGLX to read the binary file)
        new_meta_file_name = ap_meta_filename.replace('tcat', 'tcat_corrected')
        new_meta_file_path = pathlib.Path(probe_path, new_meta_file_name)
        with open(new_meta_file_path, 'w') as f:
            for key, val in ap_meta_dict.items():
                if key in ['imroTbl', 'muxTbl', 'snsChanMap', 'snsShankMap', 'snsGeomMap']:
                    f.write('~{}={}\n'.format(key, val))
                else:
                    f.write('{}={}\n'.format(key, val))

        # Read the binary data as a memory-mapped file
        ap_raw_data = readSGLX.makeMemMapRaw(ap_bin_path, ap_meta_dict)

        # Run TPrime to get artifact times aligned to probe timebase
        nidq_stream_idx = 10 # arbitrary index number
        syncperiod = config['tprime']['syncperiod']
        tostream_probe_edges_file = '{}_tcat.imec{}.ap.xd_{}_6_500.txt'.format(run_name, probe_id, int(ap_meta_dict['nSavedChans'])-1)

        command = ['TPrime',
                     '-syncperiod={}'.format(syncperiod),
                     '-tostream={}'.format(os.path.join(probe_path, tostream_probe_edges_file)),
                     '-fromstream={},{}'.format(nidq_stream_idx, os.path.join(input_dir, epoch_name, run_name + '_tcat.nidq.xa_0_0.txt')),
                     '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, epoch_name, '{}_tcat.nidq.xa_3_0.txt'.format(run_name)),
                                  os.path.join(probe_path, 'whisker_stim_times_to_imec{}.txt'.format(probe_id))),
             ]

        subprocess.run(command, shell=True, cwd=config['tprime']['tprime_path'])

        # Read the artifact times, and convert times to indices
        artifact_times = np.loadtxt(os.path.join(probe_path, 'whisker_stim_times_to_imec{}.txt'.format(probe_id)))
        fs = float(ap_meta_dict['imSampRate'])
        indices = artifact_times * fs
        indices = np.round(indices).astype(int, order='F')

        # Ensures indices are within bounds of the array
        indices = np.clip(indices, 0, ap_raw_data.shape[1] - 1)

        # Get chunks of indices with correction window
        window_ms = config['artifact_correction']['window_ms']
        window_samples = int(window_ms * fs / 1000)
        indices = np.array([np.arange(i, i + window_samples) for i in indices])

        # Create replacement values of same size as indices and number of channels
        new_values = np.zeros((ap_raw_data.shape[0], indices.shape[0], window_samples))

        # Get new values for all indices
        for i, indices_chunk in enumerate(indices):
            # Get data just before the artifact indices
            before_indices = indices_chunk - window_samples
            # Get mean of data before the artifact indices
            ch_means = np.mean(ap_raw_data[:, before_indices], axis=1)
            # Broadcast to window size
            new_values[:,i] = np.repeat(ch_means, window_samples).reshape(-1, window_samples)

        # Create a copy of the memmap file, if it does not exist
        new_ap_file_name = ap_bin_filename.replace('tcat', 'tcat_corrected')
        new_file_path = pathlib.Path(probe_path, new_ap_file_name)

        #if not os.path.exists(new_file_path):
        shutil.copyfile(ap_bin_path, new_file_path)

        # Open the memmap file in read-write mode
        data = np.memmap(new_file_path, dtype=ap_raw_data.dtype, mode='r+', shape=ap_raw_data.shape, order='F') # 'F' means column-major order (Fortran-like)

        # Debug: check what underlying data is
        debug = False
        if debug:
            fig, axs = plt.subplots(1,4, sharey=True)
            fig.suptitle('Probe {}'.format(probe_id))
            ids_to_plot_0 = indices[0]
            ids_to_plot_1 = indices[50]
            ids_to_plot_2 = indices[150]
            ids_to_plot_3 = indices[-1]
            print('Indices to plot:', ids_to_plot_1, ids_to_plot_2)
            for i in range(0,ap_raw_data.shape[0])[::10]:
                axs[0].plot(ap_raw_data[i, ids_to_plot_0], label='early')
                axs[1].plot(ap_raw_data[i, ids_to_plot_1], label='late')
                axs[2].plot(ap_raw_data[i, ids_to_plot_2], label='late')
                axs[3].plot(ap_raw_data[i, ids_to_plot_3], label='late')

            plt.show()


        # Replace data with set of all new values
        data[:, indices] = new_values

        # Write data to disk
        data.flush()
    return
