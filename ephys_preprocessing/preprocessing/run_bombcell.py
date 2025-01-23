#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_bombcell.py
@time: 4/24/2024 4:25 PM
"""

import os
import sys
import pathlib
import re
from loguru import logger
from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording
from ephys_preprocessing.utils.phylib_utils import load_phy_model

# os.environ["MATLAB_ENGINE"] = "R2021b"
import matlab.engine

def extract_ks_version(s):
    match = re.search(r'kilosort(\d+(?:\.\d+)?)', s)
    if match:
        version_str = match.group(1)
        return float(version_str) if '.' in version_str else int(version_str)
    return None


def main(input_dir, config):
    """
    Run bombcell from MATLAB on kilosort output data.
    This computes quality metrics for each cluster identified by Kilosort.
    This does not need to be run on synchronized data so it could be run after Kilosort.
    This opens Phy GUI to generate cluster_info table.
    :param input_dir:
    :param config:
    :return:
    """

    input_dir = os.path.join(input_dir, [f for f in os.listdir(input_dir) if 'catgt' in f][0])
    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Perform computations for each probe separately
    for probe_id in probe_ids:

        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[0]
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        probe_path = os.path.join(input_dir, probe_folder)

        kilosort_folders = pathlib.Path(probe_path).glob('kilosort*')
        for kilosort_folder in kilosort_folders:
            # Set paths
            ks_name = kilosort_folder.name
            kilosort_version = extract_ks_version(ks_name)

            kilosort_path = os.path.join(kilosort_folder, 'sorter_output')
            # apbin_fname = '{}_tcat_corrected.imec{}.ap.bin'.format(epoch_name, probe_id)
            apbin_fname = '{}_tcat.imec{}.ap.bin'.format(epoch_name, probe_id)
            path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)
            # meta_fname = '{}_tcat_corrected.imec{}.ap.meta'.format(epoch_name, probe_id)
            meta_fname = '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id)
            path_to_meta = os.path.join(input_dir, probe_folder, meta_fname)

            # Start MATLAB engine
            sys.path.append(config['bombcell']['matlab_path'])
            logfile_path = os.path.join(probe_path, 'run_bombcell_log.txt')
            eng = matlab.engine.start_matlab("-nodesktop -logfile " + str(logfile_path))
            eng.addpath(eng.genpath(config['bombcell']['npy_path']), nargout=0)
            eng.cd(config['bombcell']['run_bombcell_path'], nargout=0)

            logger.info('Running bombcell for IMEC probe {}.'.format(probe_id))
            eng.run_bombcell(kilosort_path, path_to_apbin, path_to_meta, kilosort_version, nargout=0)

            # Stop MATLAB engine
            eng.quit()

        # cluster_info table creation
        logger.info('Creating cluster_info table for IMEC probe {}.'.format(probe_id))
        phy_model = load_phy_model(os.path.join(kilosort_path, 'params.py'))
        phy_model.create_metrics_dataframe()
        phy_model.save_metrics_tsv(os.path.join(kilosort_path, 'cluster_info.tsv'))

    return
