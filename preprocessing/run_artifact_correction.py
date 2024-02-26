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

import numpy as np
import pathlib

from matplotlib import pyplot as plt

# Import readers
import readSGLX


def main(input_dir, config):
    """
    Run artifact correction on CatGT-processed ephys data using artifact times.
    :param input_dir:
    :param config:
    :return:
    """

    epoch_name = os.listdir(input_dir)[0]
    run_name = epoch_name.split('_')[1] + '_' + epoch_name.split('_')[2]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    n_probes = len(probe_folders)

    # Get artifact times (whisker stim times) obtained from CatGT
    path_to_artifact_times = os.path.join(input_dir, epoch_name, 'sync_event_times', 'whisker_stim_times.txt')
    artifact_times = np.loadtxt(path_to_artifact_times)

    print('Running artifact correction...')
    for probe_id in range(n_probes)[1:2]:
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

        # Convert artifact times to indices
        fs = float(ap_meta_dict['imSampRate'])
        indices = artifact_times * fs
        indices = np.round(indices).astype(int, order='F')

        # Ensures indices are within bounds of the array #useful?
        indices = np.clip(indices, 0, ap_raw_data.shape[1] - 1)

        # Get chunks of indices with correction window
        window_ms = config['window_ms']
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
        if not os.path.exists(new_file_path):
            shutil.copyfile(ap_bin_path, new_file_path)
        #shutil.copyfile(ap_bin_path, new_file_path)

        data = np.memmap(new_file_path, dtype=ap_raw_data.dtype, mode='r+', shape=ap_raw_data.shape, order='F') # 'F' means column-major order (Fortran-like)

        # Replace data with set of all new values - TEST
        data[:, indices] = new_values

        # Iterate over each index list and replace with new value
        #for idx, (indices, new_value) in enumerate(zip(indices[0:5], new_values[0:5])):
        #    print(idx, indices, np.array(indices.shape), new_value.shape)
        #    data[:, indices] = new_value.T

        # Write data to disk
        data.flush()



    return
