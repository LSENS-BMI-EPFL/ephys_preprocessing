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

    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))
    start_time = time.time()

    # Run IBL ephys-atlas formatting
    logger.info('Starting IBL ephys-atlas data formatting.')
    run_ibl_ephys_atlas_format.main(input_dir, config['anatomy'])
    logger.info('Finished IBL ephys-atlas formatting in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    return