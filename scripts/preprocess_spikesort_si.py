# Imports
import argparse
import os
import yaml
import pathlib
import time
from loguru import logger
from pathlib import Path
logger.add("log/preprocess_spikesort_{time}.log", colorize=True,
           format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
from ephys_preprocessing.preprocessing import (
    run_catgt, 
    run_sorter,
    run_bombcell,
)

@logger.catch
def main(input_dir, config_file):
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))
    start_time = time.time()

    # Get epoch number and run name
    epoch_name = os.listdir(input_dir)[0]
    # Find only directory names with "imec" in it
    n_probes = len([f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f])
    logger.info('Recording using {} probe(s).'.format(n_probes))

    # Create output folder
    # mouse_name = input_dir.split('\\')[2]
    # session_name = input_dir.split('\\')[-2]
    mouse_name = input_dir.parents[2].name
    session_name = input_dir.name
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    logger.info('Processed data will be saved to {}.'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    # Run CatGT
    logger.info('Starting CatGT.')
    run_catgt.main(input_dir, processed_dir, config['catgt'])
    logger.info('Finished CatGT in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Run Kilosort
    logger.info('Starting Kilosort.')
    run_sorter.main(processed_dir, config)
    logger.info("Finished Kilosort in {}.".format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))


    # Run quality metrics e.g. bombcell
    logger.info('Starting bombcell quality metrics.')
    run_bombcell.main(processed_dir, config)
    logger.info('Finished bombcell quality metrics in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--input', type=str, nargs='?', required=True)
    # parser.add_argument('--config', type=str, nargs='?', required=False)
    # args = parser.parse_args()

    # args = {}

    # args.input = Path('/mnt/lsens/data/PB191/Recording/Ephys/PB191_20241210_110601') #until \Ephys
    # args.config = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/ephys_preprocessing/preprocessing/preprocess_config_si.yaml')
    
    data_path = Path('/mnt/lsens/data')
    input_list = [
        'PB191/Recording/Ephys/PB191_20241210_110601',
        # 'PB192/Recording/Ephys/PB192_20241211_113347',
        # 'PB201/Recording/Ephys/PB201_20241212_192123',
        # 'PB195/Recording/Ephys/PB195_20241214_114010',
        # 'PB196/Recording/Ephys/PB196_20241217_144715',
    ]

    config = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/ephys_preprocessing/preprocessing/preprocess_config_si.yaml')

    for input in input_list:
        main(data_path / input, config)