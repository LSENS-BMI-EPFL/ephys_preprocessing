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
from pathlib import Path
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
    # logger.info('Starting LFP analysis.')
    # run_lfp_analysis.main(input_dir)
    # logger.info('Finished LFP analysis in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    exec_time_hhmmss = time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))
    logger.success(f'Finished preprocessing in {exec_time_hhmmss} for {input_dir}.')

    return


if __name__ == '__main__':
        # parser = argparse.ArgumentParser()
        # parser.add_argument('--input', type=str, nargs='?', required=True)
        # parser.add_argument('--config', type=str, nargs='?', required=False)
        # args = parser.parse_args()

        # args.input = r'M:\\analysis\\Axel_Bisi\\data\AB141\AB141_20241127_140308\Ephys\catgt_AB141_g1'
        # args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        # main(args.input, args.config)

        data_path = Path('/mnt/lsens/analysis/Jules_Lebert/data')
        input_list = [
            # 'PB191/PB191_20241210_110601/Ephys/catgt_PB191_g0',
            # 'PB192/PB192_20241211_113347/Ephys/catgt_PB192_g2',
            # 'PB201/PB201_20241212_192123/Ephys/catgt_PB201_g0',
            # 'PB195/PB195_20241214_114010/Ephys/catgt_PB195_g0',
            # 'PB196/PB196_20241217_144715/Ephys/catgt_PB196_g0',
            # 'RD076/RD076_20250214_125235/Ephys/catgt_RD076_g1',
            # 'RD077/RD077_20250219_183425/Ephys/catgt_RD077_g0',
            # 'RD077/RD077_20250221_102024/Ephys/catgt_RD077_g0',
            # 'RD072/RD072_20250305_131521/Ephys/catgt_RD072_g0',
            # 'JL002/JL002_20250507_135553/Ephys/catgt_JL002_20250507_g0',
            # 'PB193/PB193_20241218_135125/Ephys/catgt_PB193_g0',
            # 'PB194/PB194_20241218_161235/Ephys/catgt_PB194_g0',
            # 'PB197/PB197_20241216_155436/Ephys/catgt_PB197_g0',
            # 'PB198/PB198_20241213_142448/Ephys/catgt_PB198_g2',
            # 'PB200/PB200_20241216_112934/Ephys/catgt_PB200_g0',
            # TODO: Check paths below if correct
            'JL005/JL005_20250520_142542/Ephys/catgt_JL005_20250520_g0',
            'JL002/JL002_20250522_111333/Ephys/catgt_JL002_20250522_2_g0',
            'JL002/JL002_20250523_101907/Ephys/catgt_JL002_20250523_g0',
            'JL007/JL007_20250523_144849/Ephys/catgt_JL007_20250523_g0',
            'JL006/JL006_20250601_104051/Ephys/catgt_JL006_20250601_g0',
            'JL006/JL006_20250602_122916/Ephys/catgt_JL006_20250602_g0',
            'JL007/JL007_20250603_150143/Ephys/catgt_JL007_20250603_2_g0',
            'JL007/JL007_20250605_145217/Ephys/catgt_JL007_20250605_2_g0'
        ]

        config = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/scripts/preprocess_config_si.yaml')

        for input in input_list:
            main(data_path / input, config)