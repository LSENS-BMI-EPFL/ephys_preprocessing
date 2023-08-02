#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: preprocessing.py
@time: 25.04.2023 10:18
"""

# Imports
import argparse
import os
import yaml
import pathlib

import run_catgt, run_kilosort

def main(input_dir, config_file):

    with open(config_file) as file:
        config = yaml.safe_load(file)

    print('Preprocessing data from {}...'.format(input_dir))

    # Get epoch number and run name
    epoch_name = os.listdir(input_dir)[0]
    dirnames = 1  # takes first run i.e. gN_1st
    n_probes = len(next(os.walk(os.path.join(input_dir, epoch_name)))[dirnames])
    print('Recording using {} probe(s)'.format(n_probes))

    # Create output folder
    mouse_name = input_dir.split('\\')[-4]
    session_name = input_dir.split('\\')[-2]
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    print('Saving processed data to {}:'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    # Run CatGT
    run_catgt.main(input_dir, processed_dir, config['catgt'])
    print('Finished CatGT.')

    # Run Kilosort
    run_kilosort.main(processed_dir, config['kilosort'])
    print('Finished Kilosort.')

    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, nargs='?', required=True)
    parser.add_argument('--config', type=str, nargs='?', required=False)
    args = parser.parse_args()

    args.input = r'M:\data\AB077\Recording\AB077_20230531_143839\Ephys'
    args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

    main(args.input, args.config)