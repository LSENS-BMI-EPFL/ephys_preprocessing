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

import run_tprime
import run_cwaves
import run_mean_waveform_metrics
import run_lfp_analysis
from preprocessing import run_bombcell


def main(input_dir, config_file):
    """
    Run TPrime and Cwaves on preprocessed spike data.
    :param input_dir: path to CatGT preprocessed data
    :param config_file: path to config file
    :return:
    """

    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    print('Preprocessing data from {}...'.format(input_dir))

    # Run TPrime
    print('Starting Tprime.')
    run_tprime.main(input_dir, config['tprime'])
    print('Finished Tprime.')

    # Run Cwaves
    print('Starting Cwaves.')
    run_cwaves.main(input_dir, config['cwaves'])
    print('Finished Cwaves.')

    # Run mean waveform metrics
    print('Starting mean waveform metrics.')
    run_mean_waveform_metrics.main(input_dir)
    print('Finished mean waveform metrics.')

    # LFP analysis for depth estimation
    print('Starting LFP analysis.')
    run_lfp_analysis.main(input_dir)
    print('Finished LFP analysis.')

    print('Finished sync preprocessing for {}.'.format(input_dir))

    return


if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', type=str, nargs='?', required=True)
        parser.add_argument('--config', type=str, nargs='?', required=False)
        args = parser.parse_args()

        args.input = r'M:\\analysis\\Axel_Bisi\\data\AB122\AB122_20240804_134554\Ephys\catgt_AB122_g0'
        args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        main(args.input, args.config)