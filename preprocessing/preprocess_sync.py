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
    run_tprime.main(input_dir, config['tprime'])
    print('Finished Tprime.')

    # Run Cwaves
    run_cwaves.main(input_dir, config['cwaves'])
    print('Finished Cwaves.')

    print('Finished preprocessing for {}.'.format(input_dir))

    return


if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', type=str, nargs='?', required=True)
        parser.add_argument('--config', type=str, nargs='?', required=False)
        args = parser.parse_args()

        args.input = r'M:\analysis\Axel_Bisi\data\AB077\AB077_20230531_143839\Ephys\catgt_AB077_g2'
        args.config = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'

        main(args.input, args.config)