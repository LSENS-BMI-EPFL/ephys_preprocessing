#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_artifact_correction.py
@time: 2/14/2024 11:32 AM
"""


# # Imports
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import shutil
import subprocess
import numpy as np
import pathlib
from loguru import logger
from utils.ephys_utils import check_if_valid_recording

from matplotlib import pyplot as plt

# Import readers
from utils import readSGLX


def main(input_dir, config):
    """
    Run artifact correction on CatGT-processed ephys data using TPrime-aligned artifact times.
    This runs T-Prime to get artifact times aligned to probe timebase, and then replaces the artifact times with
    the mean of the data just before the artifact times.
    This reduces saturation/neuron-like/extra-filtering artifacts in the data, and is beneficial before spike sorting.
    :param input_dir:
    :param config:
    :return:
    """

    epoch_name = os.listdir(input_dir)[0]
    run_name = epoch_name[6:]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]

    for probe_folder in probe_folders:
        probe_id = int(probe_folder.split('imec')[-1])

        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[1]
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

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
        logger.info('TPrime pass to sync whisker artifact times to IMEC probe {} timebase.'.format(probe_id))
        subprocess.run(command, shell=True, cwd=config['tprime']['tprime_path'])

        # Read artifact times
        artifact_times = np.loadtxt(os.path.join(probe_path, 'whisker_stim_times_to_imec{}.txt'.format(probe_id)))
        fs = float(ap_meta_dict['imSampRate'])
        indices = np.clip((artifact_times * fs).round().astype(int), 0, ap_raw_data.shape[1] - 1)

        # Compute correction window
        window_samples = int(config['artifact_correction']['window_ms'] * fs / 1000)
        indices = np.array([np.arange(i, i + window_samples) for i in indices])

        # Optimize by replacing only affected segments
        before_indices = np.maximum(indices - window_samples, 0)

        # Compute mean of pre-artifact data (vectorized)
        ch_means = np.mean(ap_raw_data[:, before_indices], axis=2, keepdims=True)

        # Generate new values and replace in-place
        new_values = np.tile(ch_means, (1, 1, window_samples))

        # Create a copy of the memmap file, if it does not exist
        new_ap_file_name = ap_bin_filename.replace('tcat', 'tcat_corrected')
        new_ap_file_path = pathlib.Path(probe_path, new_ap_file_name)
        logger.info('Writing a corrected copy of the .ap.bin file.')
        shutil.copyfile(ap_bin_path, new_ap_file_path)  # Always overwrite

        data = np.memmap(new_ap_file_path, dtype=ap_raw_data.dtype, mode='r+', shape=ap_raw_data.shape, order='F')

        for i in range(new_values.shape[1]):  # batch-process
            idx = indices[i]  # Get index chunk
            data[:, idx] = new_values[:, i]  # vectorized batch write

        data.flush()  # Ensures all writes are committed
        del data  # Close memmap explicitly to avoid memory issues
        logger.info('Artifact correction completed for probe {}.'.format(probe_id))

    return
