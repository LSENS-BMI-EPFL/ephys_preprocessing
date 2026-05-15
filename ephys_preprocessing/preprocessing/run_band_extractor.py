#! /usr/bin/env python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_band_extractor.py
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pathlib
from loguru import logger
from neuropixel import NP2Converter


def main(input_dir):
    """
    Extract LFP band from a SpikeGLX NP2 .ap.bin file.
    :param input_dir: path to folder containing .ap.bin and .ap.meta files
    :return:
    """
    input_dir = pathlib.Path(input_dir)
    bin_file = next(input_dir.glob('*.ap.bin'))
    logger.info('AP bin file: {}'.format(bin_file.name))

    converter = NP2Converter(ap_file=bin_file, #full band file is ap.bin
                             post_check=True,
                             delete_original=False,
                             compress=False)
    logger.info('Probe version: {}'.format(converter.np_version))

    status = converter.process(overwrite=False)
    logger.info('NP2Converter status: {}'.format(status))
    return


if __name__ == '__main__':
    input_dir = sys.argv[1]
    main(input_dir=input_dir)
