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
import time
from loguru import logger
logger.add("log/preprocess_spikesort_{time}.log", colorize=True,
           format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
import run_catgt
import run_artifact_correction
import run_overstrike
import run_kilosort
import run_bombcell


@logger.catch
def main(input_dir, config_file):
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))

    # Get epoch number and run name
    epoch_name = os.listdir(input_dir)[0]
    # Find only directory names with "imec" in it
    n_probes = len([f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f])
    logger.info('Recording using {} probe(s).'.format(n_probes))

    # Create output folder
    mouse_name = input_dir.split('\\')[2]
    session_name = input_dir.split('\\')[-2]
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    logger.info('Processed data will be saved to {}.'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    # Run CatGT
    logger.info('Starting CatGT.')
    run_catgt.main(input_dir, processed_dir, config['catgt'])
    logger.info('Finished CatGT.')

    # Run TPrime a first time to sync whisker artifact times
    logger.info('Starting artifact correction.')
    run_artifact_correction.main(processed_dir, config)
    logger.info('Finished artifact correction.')

    # Optionally, run OverStrike
    perform_overstrike = False
    if perform_overstrike:

        # List of time spans to zero out in recording in secs, relative to start of recording
        if mouse_name == 'AB105':
            timespans_list = [(0,24),(1800,1933),(2389,3161)]
        elif mouse_name == 'AB107':
            timespans_list = [(1247, 1930)]
        #elif mouse_name == 'AB129':
        #    timespans_list = [(5580, 1E9)]
        elif mouse_name == 'AB142':
            timespans_list = [(3474, 3892)]
        else:
            timespans_list = [(), ]


        # Run overstrike on all probes
        logger.info('Starting OverStrike.')
        run_overstrike.main(processed_dir, config['overstrike'], timespans_list=timespans_list)
        logger.info('Finished OverStrike.')

    # Run Kilosort
    logger.info('Starting Kilosort.')
    #run_kilosort.main(processed_dir, config)
    logger.info("Finished Kilosort.")

    # Run quality metrics e.g. bombcell
    logger.info('Starting bombcell quality metrics.')
    #run_bombcell.main(processed_dir, config)
    logger.info('Finished bombcell quality metrics.')

    catgt_epoch_name = os.listdir(processed_dir)[0]
    logger.success(f'Finished preprocessing & spike sorting in: {os.path.join(processed_dir, catgt_epoch_name)}. You '
                   f'can now visually check spike sorting results using Phy, then use this path as input to the '
                   f'script preprocessing_sync.py.')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, nargs='?', required=True)
    parser.add_argument('--config', type=str, nargs='?', required=False)
    args = parser.parse_args()

    args.input = r'M:\data\AB143\Recording\AB143_20241126_115737\Ephys' #until \Ephys
    args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

    main(args.input, args.config)
