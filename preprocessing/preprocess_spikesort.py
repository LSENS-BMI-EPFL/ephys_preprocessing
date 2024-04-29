#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: preprocess_spikesort.sbatch.py
@time: 25.04.2023 10:18
"""

# Imports
import argparse
import os
import yaml
import pathlib

# Import submodules
import run_catgt
import run_tprime
import run_artifact_correction
import run_overstrike
import run_kilosort
from preprocessing import run_bombcell


def main(input_dir, config_file):
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    print('Preprocessing data from {}...'.format(input_dir))

    # Get epoch number and run name
    epoch_name = os.listdir(input_dir)[0]
    # Find only directory names with "imec" in it
    n_probes = len([f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f])
    print('Recording using {} probe(s)'.format(n_probes))

    # Create output folder
    mouse_name = input_dir.split('\\')[2]
    session_name = input_dir.split('\\')[-2]
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    print('Saving processed data to {}:'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    # Run CatGT
    #run_catgt.main(input_dir, processed_dir, config['catgt'])
    print('Finished CatGT.')

    # Run TPrime a first time to sync whisker artifact times
    #run_artifact_correction.main(processed_dir, config)
    print('Finished artifact correction.')

    # Optionally, run OverStrike
    perform_overstrike=False
    if perform_overstrike:

        # List of time spans to zero out in recording in secs, relative to start of recording
        if mouse_name == 'AB105':
            timespans_list = [(0,24),(1800,1933),(2389,3161)]
        elif mouse_name == 'AB107':
            timespans_list = [(1247, 1930)]
        else:
            timespans_list = [(), ]

        # Run overstrike on all probes
        run_overstrike.main(processed_dir, config['overstrike'], timespans_list=timespans_list)
        print('Finished OverStrike.')

    # Run Kilosort
    #run_kilosort.main(processed_dir, config)
    print('Finished Kilosort.')

    # Run quality metrics e.g. bombcell
    run_bombcell.main(processed_dir, config)
    print('Finished bombcell quality metrics.')


    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, nargs='?', required=True)
    parser.add_argument('--config', type=str, nargs='?', required=False)
    args = parser.parse_args()

    #args.input = r'M:\analysis\Axel_Bisi\data\AB085\AB085_20231005_152636\Ephys' #until \Ephys
    args.input = r'M:\data\AB105\Recording\AB105_20240314_115206\Ephys'
    args.input = r'M:\data\AB102\Recording\AB102_20240309_114107\Ephys'
    args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

    main(args.input, args.config)
