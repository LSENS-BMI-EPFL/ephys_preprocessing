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
    mouse_name = input_dir.split('\\')[4]
    session_name = input_dir.split('\\')[-2]
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    print('Saving processed data to {}:'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)


    # Run CatGT
    #run_catgt.main(input_dir, processed_dir, config['catgt'])
    print('Finished CatGT.')

    # Run TPrime a first time to sync whisker artifact times
    #run_tprime.main(processed_dir, config['tprime'], pre_ks=True)
    run_artifact_correction.main(processed_dir, config['artifact_correction'])
    print('Finished TPrime and artifact correction.')

    # Optionally, run OverStrike
    perform_overstrike=False
    if perform_overstrike:

        # List of time spans to zero out in recording
        timespans_list = [(),] # in secs, relative to start of recording

        # Run overstrike on all probes
        run_overstrike.main(processed_dir, config['overstrike'],
                            timespans_list=timespans_list)


        print('Finished OverStrike.')

    # Run Kilosort
    #run_kilosort.main(processed_dir, config['kilosort'])
    #print('Finished Kilosort.')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, nargs='?', required=True)
    parser.add_argument('--config', type=str, nargs='?', required=False)
    args = parser.parse_args()

    args.input = r'M:\analysis\Axel_Bisi\data\AB092\AB092_20231205_140109\Ephys' #until \Ephys
    args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

    main(args.input, args.config)
