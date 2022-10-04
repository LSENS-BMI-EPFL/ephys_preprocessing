#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: do_kilosort.py
@time: 16/05/2022 13:59
"""

# Imports
# import matlab.engine #must be imported first on some systems
import sys
import json
import pprint
import tkinter.filedialog as fdialog
import os
from pathlib import Path

# Modules
sys.path.append('C:\\Users\\bisi\\Github\\ecephys_spike_sorting\\ecephys_spike_sorting\\common')
from SGLXMetaToCoords import readMeta

# Load config
print('Loading Kilosort config')
with open('configs/kilosort_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

# Select preprocessed ephys data
input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory', initialdir=config['data_path'])
input_dir = os.path.join(input_dir_mouse, 'Recording/Ephys')

epoch_name = os.listdir(input_dir)[0]
run_name = os.listdir(input_dir)[0][0:-3]
run_name = run_name.replace('catgt_', '')

# Get run info and number of probes/channels saved
dirnames = 1
n_probes = len(next(os.walk(os.path.join(input_dir, epoch_name)))[dirnames])
print('Recorded using {} probes'.format(n_probes))
n_saved_ch_probes = dict.fromkeys(range(n_probes))

for probe_id in range(n_probes):
    probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)

    metafile_name = '{}_tcat.imec{}.ap.meta'.format(epoch_name.replace('catgt_', ''), probe_id)
    metafile_path = os.path.join(input_dir, epoch_name, probe_folder, metafile_name)

    meta_dict = readMeta(Path(metafile_path))

    # Add number of saved channels
    print('-- IMEC probe', probe_id, ' #channels saved:', meta_dict['nSavedChans'])
    n_saved_ch_probes[probe_id] = meta_dict['nSavedChans']

    # Create output folder for each probe
    ks_v1, _, ks_v2 = str(config['ks_version']).partition('.')
    ks_output_folder = os.path.join(input_dir, epoch_name, probe_folder, 'ks{}'.format(ks_v1 + ks_v2))
    print(ks_output_folder)
    Path(ks_output_folder).mkdir(parents=True, exist_ok=True)  # create output dir

# RUN KILOSORT FOR EACH PROBE
# def call_kilosort():
#
#    eng = matlab.engine.start_matlab()
#    eng.addpath(eng.genpath(KS_dir))
#    eng.addpath(eng.genpath(NPY_dir))
#    eng.addpath(home_dir)

for probe_id in range(n_probes):
    print('- Running Kilosort for IMEC probe', probe_id)
    ### FILL IN HERE
