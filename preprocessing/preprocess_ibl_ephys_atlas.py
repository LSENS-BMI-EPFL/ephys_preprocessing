#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: preprocess_ibl_ephys_atlas.py
@time: 8/17/2025 12:24 PM
"""

# Imports
import argparse
import yaml
import time
import subprocess
from loguru import logger
logger.add("log/preprocess_ibl_ephys_atllas_{time}.log", colorize=True,
              format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
import run_ibl_ephys_atlas_format

@logger.catch
def main(input_dir, config_file):
    """
    Run IBL formatting for the ephys-atlas alignment GUI.
    :param input_dir: path to CatGT preprocessed data
    :param config_file: path to config file
    :return:
    """
    subprocess.run('conda activate iblenv', shell=True)


    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))
    start_time = time.time()

    # Run IBL ephys-atlas formatting
    logger.info('Starting IBL ephys-atlas data formatting.')
    run_ibl_ephys_atlas_format.main(input_dir, config)
    logger.info('Finished IBL ephys-atlas formatting in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Open GUI
    logger.info('Opening ephys-atlas GUI.')
    path_to_gui = config['anatomy']['gui_path']
    command = f'python {path_to_gui} -o True'
    subprocess.run(command, shell=True)

    return

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', type=str, nargs='?', required=True)
        parser.add_argument('--config', type=str, nargs='?', required=False)
        args = parser.parse_args()

        args.input = r'M:\analysis\Axel_Bisi\data\AB164\AB164_20250422_115457\Ephys\catgt_AB164_g0'
        args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        main(args.input, args.config)