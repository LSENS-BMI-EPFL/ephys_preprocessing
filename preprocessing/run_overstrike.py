#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_overstrike.py
@time: 10/25/2023 12:29 PM
"""

# Imports
import os
import subprocess
from utils.ephys_utils import flatten_list
import webbrowser
import yaml
from loguru import logger


def main(input_dir, config, timespans_list):
    """
    Run OverStrike on raw ephys data and save to output directory.
    :param input_dir: path to raw ephys data
    :param output_dir: path to output directory
    :param config: config dict
    :param timespans_list: list of timespans to zero-out
    :return:
    """

    epoch_name = os.listdir(input_dir)[0]
    probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    n_probes = len(probe_folders)

    # Check timespans_list is a list of tuples
    try:
        assert isinstance(timespans_list, list)
        assert all(isinstance(timespan, tuple) for timespan in timespans_list)

    except AssertionError:
        raise TypeError('Timespans_list must be a list of tuples. Skipping OverStrike...')
        logger.error('Timespans_list must be a list of tuples. Skipping OverStrike.')
        return

    # Check not empty
    try:
        assert len(timespans_list) > 0
    except AssertionError:
        raise ValueError('Timespans_list cannot be empty. Skipping OverStrike...')
        logger.error('Timespans_list cannot be empty. Skipping OverStrike.')
        return

    logger.info('Striking timespans: {}'.format(timespans_list))
    for probe_id in range(n_probes):
        if probe_id != 4:
            print('Skipping OverStrike of probe {}'.format(probe_id))
            continue
        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)
        ap_bin_file_name = [f for f in os.listdir(probe_path) if 'ap.bin' in f and 'corrected' in f][0]
        ap_bin_path = os.path.join(probe_path, ap_bin_file_name)

        # Iterate over timespans to zero-out
        for timespan in timespans_list:

            # Write OverStrike command line
            command = ['OverStrike',
                       '-file={}'.format(ap_bin_path),
                       '-chans=0:300',
                       '-secs={},{}'.format(timespan[0], timespan[1])
                       ]
            logger.info('OverStrike command line will run: {}'.format(list(flatten_list(command))))

            # Run OverStrike
            subprocess.run(list(flatten_list(command)), shell=True, cwd=config['overstrike_path'])

            logger.info('Opening OverStrike log file at: {}'.format(os.path.join(config['overstrike_path'], 'OverStrike.log')))
            webbrowser.open(os.path.join(config['overstrike_path'], 'OverStrike.log'))

    # Save overstrike information
    overstrike_info = {'timespans_list': timespans_list}
    with open(os.path.join(input_dir, 'overstrike_info.yaml'), 'w') as f:
        yaml.dump(overstrike_info, f)

    return

