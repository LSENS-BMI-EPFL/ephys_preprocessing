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
import run_tprime
import run_cwaves
import run_mean_waveform_metrics
import run_lfp_analysis


@logger.catch
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
    #run_tprime.main(input_dir, config['tprime'])
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
    run_lfp_analysis.main(input_dir, config)
    logger.info('Finished LFP analysis in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    exec_time_hhmmss = time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))
    logger.success(f'Finished preprocessing in {exec_time_hhmmss} for {input_dir}.')

    return


if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', type=str, nargs='?', required=True)
        parser.add_argument('--config', type=str, nargs='?', required=False)
        args = parser.parse_args()

        args.input = r'M:\analysis\Axel_Bisi\data\PB191\PB191_20241210_110601\Ephys\catgt_PB191_g0'
        #M:\analysis\Axel_Bisi\data\PB195\PB195_20241214_114010\Ephys\catgt_PB195_g0
        #args.input = r'M:\analysis\Axel_Bisi\data\PB191\PB191_20241210_110601\Ephys\catgt_PB191_g0'
        args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        main(args.input, args.config)