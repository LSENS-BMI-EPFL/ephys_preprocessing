#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_kilosort.py
@time: 8/2/2023 5:31 PM
"""
import os
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
    print('Data to sort: ', probe_folders)
    n_probes = len(probe_folders)

    print('Running Kilosort...')
    for probe_id in range(n_probes):
        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)

        # Start MATLAB engine
        eng = matlab.engine.start_matlab()
        eng.addpath(eng.genpath(r'C:\Users\bisi\Github\npy-matlab'), nargout=0)
        eng.cd(config['kilosort_path'], nargout=0)

        # Run Kilosort for current probe
        print('- Running Kilosort for IMEC probe', probe_id)
        #eng.run_main_kilosort(probe_path, config['temp_data_path'], nargout=0)

        # Stop MATLAB engin
        eng.quit()

    return
