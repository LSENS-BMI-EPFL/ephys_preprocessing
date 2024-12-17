#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: preprocess_sync.py
@time: 8/2/2023 6:28 PM
"""

# Imports
import argparse
import yaml
import time
from loguru import logger
logger.add("log/preprocess_sync_{time}.log", colorize=True,
              format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
from ephys_preprocessing.preprocessing import (
     run_tprime,
     run_cwaves,
     run_mean_waveform_metrics,
     run_lfp_analysis,
)


def main(input_dir, config_file):
    """
    Run TPrime and Cwaves on preprocessed spike data.
    :param input_dir: path to CatGT preprocessed data
    :param config_file: path to config file
    :return:
    """

    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))
    start_time = time.time()

    # Run TPrime
    logger.info('Starting Tprime.')
    run_tprime.main(input_dir, config['tprime'])
    logger.info('Finished Tprime in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Run Cwaves
    logger.info('Starting Cwaves.')
    run_cwaves.main(input_dir, config['cwaves'])
    logger.info('Finished Cwaves in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Run mean waveform metrics
    logger.info('Starting mean waveform metrics.')
    run_mean_waveform_metrics.main(input_dir)
    logger.info('Finished mean waveform metrics in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # LFP analysis for depth estimation
    logger.info('Starting LFP analysis.')
    run_lfp_analysis.main(input_dir)
    logger.info('Finished LFP analysis in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    exec_time_hhmmss = time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))
    logger.success(f'Finished preprocessing in {exec_time_hhmmss} for {input_dir}.')

    return


if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', type=str, nargs='?', required=True)
        parser.add_argument('--config', type=str, nargs='?', required=False)
        args = parser.parse_args()

        args.input = r'M:\\analysis\\Axel_Bisi\\data\AB141\AB141_20241127_140308\Ephys\catgt_AB141_g1'
        args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        main(args.input, args.config)