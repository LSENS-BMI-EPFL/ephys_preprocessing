#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: run_ibl_ephys_atlas_format.py
@time: 8/17/2025 12:10 PM
"""


# Imports
import os
import shutil
import numpy as np
import json
from loguru import logger
from pathlib import Path

from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording

from atlaselectrophysiology.extract_files import extract_data
from iblatlas.atlas import AllenAtlas

# Monkey-patch phylib to use shutil.copyfile instead of shutil.copy
# to avoid PermissionError when overwriting files on the server
import phylib.io.alf as _alf_module

_orig_create_if_possible = _alf_module._create_if_possible


def _copy_if_possible_patched(path, new_path, force=False):
    if not _orig_create_if_possible(path, new_path, force=force):
        return False
    shutil.copyfile(str(path), str(new_path))
    return True


_alf_module._copy_if_possible = _copy_if_possible_patched

def main(input_dir, config):
    """
    Run IBl code to format ephys data for the ephys-atlas alignement GUI.
    :param input_dir:  path to preprocessed data
    :param config:  config dict
    :return:
    """
    day_index = 0       # day of recordings, iterate
    input_dir = Path(input_dir)
    catgt_epoch_name = input_dir.name
    session_date = input_dir.parents[1].name
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    mouse_id = epoch_name.split('_')[0]
    anat_data_folder = os.path.join(config['anatomy']['anat_data_path'], mouse_id, 'fused',
                                    'registered', 'segmentation', 'atlas_space', 'tracks')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = sorted([f[-1] for f in probe_folders])

    output_path = config['output_path']

    # Perform computations for each probe separately
    for probe_id in probe_ids:

        if not check_if_valid_recording(config, mouse_id, probe_id, day_id=day_index):
            continue

        # Format electrophysiology data
        # -----------------------------
        logger.info(f'- IMEC {probe_id} ephys spike-sorting and output data...')
        # Path to Kilosort output
        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        ephys_path = input_dir / probe_folder

        kilosort_folders = list(ephys_path.glob('kilosort*'))
        if len(kilosort_folders) == 0:
            logger.error(f"Kilosort output not found at {kilosort_folders}. Please run Kilosort first.")
        if len(kilosort_folders) > 1:
            logger.warning(f"Not implemented: multiple Kilosort outputs found at {kilosort_folders}. Using {kilosort_folders[0]}.")

        # Path to raw ephys data
        if not ephys_path.exists():
            logger.error(f"Ephys data not found at {ephys_path}. Please check the path.")

        # Save path
        if mouse_id.startswith('AB'):
            out_path = Path(input_dir) / probe_folder / 'ibl_format'
        elif mouse_id.startswith('MH'):
            out_path = Path(output_path) / mouse_id / session_date / 'Ephys' / catgt_epoch_name / probe_folder / 'ibl_format'
        else:
            out_path = ephys_path / 'ibl_format'

        xyz_picks_path = out_path / 'xyz_picks.json'
        overwrite = config['anatomy']['overwrite']

        if xyz_picks_path.exists() and not overwrite:
            logger.info(f'ibl_format folder already exists at {out_path}. Skipping ibl format conversion.')
            continue
        # If json does not exist or overwrite is True, but output dir exists, delete it before proceeding
        # if out_path.exists():
        #     logger.warning(f'Removing existing out_path directory at {out_path} before re-creating ibl format outputs.')
        #     import shutil
        #     shutil.rmtree(out_path)

        ks_path = kilosort_folders[0] / 'sorter_output'
        if not out_path.exists():
            extract_data(ks_path, ephys_path, out_path)

        # Format anatomical data from Brainreg track tracing
        # --------------------------------------------------
        logger.info(f'- IMEC {probe_id} probe track tracing data...')
        if not os.path.exists(anat_data_folder):
            logger.warning(f'Probe track tracing folder not found at {anat_data_folder}. Skipping data formatting.')
            continue

        atlas = AllenAtlas(res_um=25) # bregma estimate done in 25 micron resolution

        

        # Set path to brainreg-segment output file - day 0 vs. expert
        if mouse_id.startswith('AB'):
            brainreg_path = Path(anat_data_folder, f'imec{probe_id}.npy')
        elif mouse_id.startswith('MH'):
            if int(mouse_id[-3:]) < 20:
                brainreg_path = Path(anat_data_folder, f'imec{probe_id}.npy')
            else:
                brainreg_path = Path(anat_data_folder, session_date, f'imec{probe_id}.npy')

        # Flexible path for combining dataset with single and multi-day recordings, assumes no imec*.npy files in anat_data_folder if multi day
        else:
            # Look for imec*.npy files in anat_data_folder; if found, use them, else look in anat_data_folder / session_date
            imec_npy_files = list(Path(anat_data_folder).glob(f'imec{probe_id}.npy'))
            if len(imec_npy_files) > 0:
                brainreg_path = imec_npy_files[0]
            else:
                brainreg_path = Path(anat_data_folder, session_date, f'imec{probe_id}.npy')

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
        logger.info(f'Saved xyz_picks.json to {Path(out_path, "xyz_picks.json")}.')

    return

if __name__ == '__main__':
    import yaml
    # input_dir = Path('/mnt/lsens-data/JL006/Recording/JL006_20250601_104051/Ephys/JL006_20250601_g0')
    input_dir = Path('/mnt/lsens-analysis/Jules_Lebert/data_spikesorted/JL014/JL014_20251207_140401/Ephys/catgt_JL014_20251207_g0/')
    config_file = Path('/home/lebert/code/ephys_preprocessing/config/preprocess_config_si.yaml')
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    main(input_dir, config)