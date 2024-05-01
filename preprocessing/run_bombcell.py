#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_bombcell.py
@time: 4/24/2024 4:25 PM
"""

import os
import subprocess
import sys
import pyautogui
import time

from utils.ephys_utils import check_if_valid_recording

os.environ["MATLAB_ENGINE"] = "R2021b"
import matlab.engine

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

        # Set paths
        kilosort_folder = 'kilosort2'
        kilosort_version = 2

        kilosort_path = os.path.join(input_dir, probe_folder, kilosort_folder)
        apbin_fname = '{}_tcat_corrected.imec{}.ap.bin'.format(epoch_name, probe_id)
        path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)
        meta_fname = '{}_tcat_corrected.imec{}.ap.meta'.format(epoch_name, probe_id)
        path_to_meta = os.path.join(input_dir, probe_folder, meta_fname)

        # Start MATLAB engine
        sys.path.append(config['bombcell']['matlab_path'])
        logfile_path = os.path.join(probe_path, 'run_bombcell_log.txt')
        eng = matlab.engine.start_matlab("-logfile " + str(logfile_path))
        eng.addpath(eng.genpath(r'C:\Users\bisi\Github\npy-matlab'), nargout=0)
        eng.cd(config['bombcell']['bombcell_path'], nargout=0)

        print('- Running bombcell for IMEC probe', probe_id)
        eng.run_bombcell(kilosort_path, path_to_apbin, path_to_meta, kilosort_version, nargout=0)

        # Stop MATLAB engine
        eng.quit()

        # Execute Phy to generate cluster_info table
        command = 'conda activate phy2 && phy template-gui params.py && conda deactivate'
        process = subprocess.Popen(command,  shell=True, cwd=os.path.join(probe_path, 'kilosort2'))
        pyautogui.FAILSAFE = False # disable mouse moving fail-safe
        time.sleep(15) # wait for GUI to load
        pyautogui.hotkey('ctrl', 's') # simulate CTRL+S key press to save the file
        time.sleep(15) # wait for saving to complete
        pyautogui.hotkey('ctrl', 'q') # close GUI
        time.sleep(15)
        process.terminate() # terminate process

    return
