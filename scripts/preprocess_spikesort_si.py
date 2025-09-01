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
    run_overstrike,
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
    epoch_name = [dir for dir in os.listdir(input_dir) if not dir.startswith('.')][0]
    # Find only directory names with "imec" in it
    n_probes = len([f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f])
    logger.info('Recording using {} probe(s).'.format(n_probes))

    # Create output folder
    # mouse_name = input_dir.split('\\')[2]
    # session_name = input_dir.split('\\')[-2]
    mouse_name = input_dir.parents[2].name
    session_name = input_dir.name
    if session_name == 'Ephys':
        session_name = input_dir.parents[0].name
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    logger.info('Processed data will be saved to {}.'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    # Run CatGT
    logger.info('Starting CatGT.')
    run_catgt.main(input_dir, processed_dir, config['catgt'])
    logger.info('Finished CatGT in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Optionally, run OverStrike
    timespans_list = None
    if mouse_name == 'PB191':
        timespans_list = [(2350, 2373), (2724, 2778)]

    if timespans_list:
        logger.info('Starting OverStrike.')
        run_overstrike.main(processed_dir, config['overstrike'], timespans_list=timespans_list)
        logger.info('Finished OverStrike in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

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
        # 'PB191/Recording/Ephys/PB191_20241210_110601',
        # 'PB192/Recording/Ephys/PB192_20241211_113347',
        # 'PB201/Recording/Ephys/PB201_20241212_192123',
        # 'PB193/Recording/Ephys/PB193_20241218_135125',
        # 'PB194/Recording/Ephys/PB194_20241218_161235',
        # 'PB195/Recording/Ephys/PB195_20241214_114010',
        # 'PB196/Recording/Ephys/PB196_20241217_144715',
        # 'PB197/Recording/Ephys/PB197_20241216_155436',
        # 'PB198/Recording/Ephys/PB198_20241213_142448',
        # 'PB200/Recording/Ephys/PB200_20241216_112934',
        # 'RD076/Recording/RD076_20250214_125235/Ephys'
        # 'RD077/Recording/RD077_20250219_183425/Ephys',
        # 'RD077/Recording/RD077_20250221_102024/Ephys',
        # 'RD072/Recording/RD072_20250305_131521/Ephys',
        # 'JL002/Recording/JL002_20250507_135553/Ephys',
        # 'JL005/Recording/JL005_20250520_142542/Ephys/',
        # 'JL002/Recording/JL002_20250522_111333/Ephys/',
        # 'JL002/Recording/JL002_20250523_101907/Ephys/',
        # 'JL007/Recording/JL007_20250523_144849/Ephys/',
        # 'JL006/Recording/JL006_20250601_104051/Ephys/',
        # 'JL006/Recording/JL006_20250602_122916/Ephys/',
        'JL007/Recording/JL007_20250603_150143/Ephys/',
        'JL007/Recording/JL007_20250605_145217/Ephys/'
    ]

    config = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/scripts/preprocess_config_si.yaml')

    for input in input_list:
        main(data_path / input, config)