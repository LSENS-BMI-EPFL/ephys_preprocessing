#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_kilosort.py
@time: 8/2/2023 5:31 PM
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pathlib
from loguru import logger

from ephys_preprocessing.utils import dredge_utils
from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording


def main(input_dir, config):
    """
    Run Kilosort from MATLAB on preprocessed data.
    :param input_dir:  path to preprocessed data
    :param config:  config dict
    :return:
    """

    epoch_name = [f for f in os.listdir(input_dir) if '_g' in f][0]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    logger.info('Data to process: {}'.format(sorted(probe_folders)))
    n_probes = len(probe_folders)

    for probe_id in range(n_probes):

        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[1]
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)

        # Create output folder

        bin_file_name = [f for f in os.listdir(probe_path) if 'ap.bin' in f][0]
        bin_path = pathlib.Path(os.path.join(probe_path, bin_file_name))

        # Run DREDge pipeline
        logger.info('Running DREDge pipeline on probe {}.'.format(probe_id))
        preset=config['preset']
        out_path = pathlib.Path(os.path.join(probe_path, preset))
        dredge_utils.run(bin_file=bin_path, output_folder=out_path, preset=preset, use_lfp=False, overwrite=True)
        #out_path = pathlib.Path(os.path.join(probe_path, preset'))
        #dredge_utils.run(bin_file=bin_path, output_folder=out_path, preset='dredge_fast', use_lfp=False, overwrite=True)

    return


if __name__ == '__main__':
    main()
