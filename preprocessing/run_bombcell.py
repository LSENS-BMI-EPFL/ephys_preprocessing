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
os.environ["MATLAB_ENGINE"] = "R2021b"
import matlab.engine

def main(input_dir, config):
    """
    Run bombcell from MATLAB on kilosort output data.
    :param input_dir:
    :param config:
    :return:
    """

    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Perform computations for each probe separately
    for probe_id in probe_ids:

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        probe_path = os.path.join(input_dir, probe_folder)

        # Set paths
        kilosort_folder = [f for f in os.listdir(probe_path) if 'kilosort' in f][0] # get kilosort folder available
        kilosort_version = int(kilosort_folder[-1])
        kilosort_path = os.path.join(input_dir, probe_folder, kilosort_folder)
        apbin_fname = '{}_tcat_corrected.imec{}.ap.bin'.format(epoch_name, probe_id)
        path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)
        meta_fname = '{}_tcat_corrected.imec{}.ap.meta'.format(epoch_name, probe_id)
        path_to_meta = os.path.join(input_dir, probe_folder, meta_fname)

        # Start MATLAB engine
        sys.path.append(config['matlab_path'])
        logfile_path = os.path.join(probe_path, 'run_bombcell_log.txt')
        eng = matlab.engine.start_matlab("-logfile " + str(logfile_path))
        eng.addpath(eng.genpath(r'C:\Users\bisi\Github\npy-matlab'), nargout=0)
        eng.cd(config['bombcell_path'], nargout=0)

        print('- Running bombcell for IMEC probe', probe_id)
        eng.run_bombcell(kilosort_path, path_to_apbin, path_to_meta, kilosort_version, nargout=0)

        # Stop MATLAB engine
        eng.quit()

    return
