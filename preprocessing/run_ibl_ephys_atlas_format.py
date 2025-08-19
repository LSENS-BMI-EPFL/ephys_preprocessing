#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: run_ibl_ephys_atlas_format.py
@time: 8/17/2025 12:10 PM
"""


# Imports
import os
import numpy as np
import json
from loguru import logger
from pathlib import Path

import sys
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
from utils.ephys_utils import check_if_valid_recording
from atlaselectrophysiology.extract_files import extract_data
from iblatlas.atlas import AllenAtlas

def main(input_dir, config):
    """
    Run IBl code to format ephys data for the ephys-atlas alignement GUI.
    :param input_dir:  path to preprocessed data
    :param config:  config dict
    :return:
    """

    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    mouse_id = epoch_name.split('_')[0]
    anat_data_folder = os.path.join(config['anatomy']['anat_data_path'], mouse_id, 'fused',
                                    'registered', 'segmentation', 'atlas_space', 'tracks')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = sorted([f[-1] for f in probe_folders])

    # Perform computations for each probe separately
    for probe_id in probe_ids:


        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        # Format electrophysiology data
        # -----------------------------
        logger.info('- Ephys spike-sorting and output data...')
        # Path to Kilosort output
        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        ks_path = Path(os.path.join(input_dir, probe_folder, 'kilosort2'))
        if not ks_path.exists():
            logger.error(f"Kilosort output not found at {ks_path}. Please run Kilosort first.")

        # Path to raw ephys data
        ephys_path = Path(input_dir, probe_folder)
        if not ephys_path.exists():
            logger.error(f"Ephys data not found at {ephys_path}. Please check the path.")

        # Save path
        out_path = Path(os.path.join(input_dir, probe_folder, 'ibl_format'))
        if not out_path.exists():
            out_path.mkdir(parents=True, exist_ok=True)

        extract_data(ks_path, ephys_path, out_path)

        # Format anatomical data from Brainreg track tracing
        # --------------------------------------------------
        logger.info('- Probe track tracing data...')
        if not os.path.exists(anat_data_folder):
            logger.warning(f'Probe track tracing folder not found at {anat_data_folder}. Skipping data formatting.')
            continue

        atlas = AllenAtlas(res_um=25) # bregma estimate done in 25 micron resolution

        brainreg_path = Path(anat_data_folder, f'imec{probe_id}.npy')
        if not brainreg_path.exists():
            logger.warning(f'Brainreg-segment track file not found for IMEC {probe_id}. Check folder.')
            continue

        # Load in coordinates of track in CCF space (order - apdvml, origin - top, left, front voxel
        xyz_apdvml = np.load(brainreg_path)

        # Convert to IBL space (order - mlapdv, origin - bregma)
        xyz_mlapdv = atlas.ccf2xyz(xyz_apdvml, ccf_order='apdvml') * 1e6
        xyz_picks = {'xyz_picks': xyz_mlapdv.tolist()}

        # Path to save the data (same folder as where you have the ephys data)
        with open(Path(out_path, 'xyz_picks.json'), "w") as f:
            json.dump(xyz_picks, f, indent=2)

    return

