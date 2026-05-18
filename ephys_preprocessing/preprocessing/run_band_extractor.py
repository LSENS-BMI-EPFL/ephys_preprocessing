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
import re
from loguru import logger
from neuropixel import NP2Converter


from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording, get_probe_version, get_exp_datetime


def main(input_dir, config):
    """
    Extract LFP band from a SpikeGLX NP2 .ap.bin file.
    :param input_dir: path to folder containing .ap.bin and .ap.meta files
    :return:
    """

    input_dir = pathlib.Path(input_dir)
    date = get_exp_datetime(input_dir)
    epoch_name = next(d for d in input_dir.iterdir() if d.is_dir()).name
    #probe_folders = sorted((input_dir / epoch_name).glob('*imec[0-9]*'))

    probe_folders = sorted(
        p for p in (input_dir / epoch_name).iterdir()
        if p.is_dir() and re.fullmatch(r'.*imec\d+', p.name)            # this looks exactly for the original probe folders e.g. imec0, imec1
    )                                                                           # will exclude imec1a folders that are already split
    logger.info('Found {} probe(s): {}'.format(len(probe_folders), [p.name for p in probe_folders]))

    mouse_id = epoch_name.split('_')[0]

    for probe_folder in probe_folders:
        probe_id = int(probe_folder.name.split('imec')[-1])
        if not check_if_valid_recording(config, mouse_id, probe_id, date):
            continue
        if get_probe_version(config, mouse_id, probe_id, date) == 1:
            logger.info(f"Probe {probe_id} for mouse {mouse_id} is version NP1, skipping band/shank data extraction.")
            continue
        elif get_probe_version(config, mouse_id, probe_id, date) == 2:
            logger.info(f"Probe {probe_id} for mouse {mouse_id} is version NP2, proceeding with LFP band/shank data extraction.")

        bin_file = next(probe_folder.glob('*.imec*.ap.bin'))
        logger.info('Processing: {}'.format(bin_file.name))
        converter = NP2Converter(ap_file=bin_file,
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
