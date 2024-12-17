#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_kilosort.py
@time: 8/2/2023 5:31 PM
"""
import os
import sys
import pathlib
from loguru import logger
from ephys_preprocessing.utils import readSGLX
from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording

os.environ["MATLAB_ENGINE"] = "R2021b"
import matlab.engine


def main(input_dir, config):
    """
    Run Kilosort from MATLAB on preprocessed data.
    :param input_dir:  path to preprocessed data
    :param config:  config dict
    :return:
    """

    epoch_name = os.listdir(input_dir)[0]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    logger.info('Data to spike-sort: {}'.format(probe_folders))
    n_probes = len(probe_folders)

    for probe_id in range(n_probes):

        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[1]
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)

        # Create output folder
        pathlib.Path(os.path.join(probe_path, 'kilosort2')).mkdir(parents=True, exist_ok=True)

        meta_file_name = [f for f in os.listdir(probe_path) if 'ap.meta' in f][0]
        ap_meta_config = readSGLX.readMeta(pathlib.Path(probe_path, meta_file_name))
        fs = float(ap_meta_config['imSampRate'])

        # Start MATLAB engine
        sys.path.append(config['kilosort']['matlab_path'])
        logfile_path = os.path.join(probe_path, 'preprocess_spikesort_log.txt')
        eng = matlab.engine.start_matlab("-logfile " + str(logfile_path))
        eng.cd(config['kilosort']['kilosort_path'], nargout=0)

        logger.info('Running Kilosort for IMEC probe {}.'.format(probe_id))
        eng.run_main_kilosort(probe_path, fs, config['kilosort']['temp_data_path'], nargout=0)

        # Stop MATLAB engine
        eng.quit()

    return
