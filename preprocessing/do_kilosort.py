#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: do_kilosort.py
@time: 16/05/2022 13:59
"""

# Imports
# import matlab.engine #must be imported first on some systems
import json
import pprint
from tkinter.filedialog import askdirectory
import os
import pathlib
import readSGLX
import matlab.engine

# Load config
print('Loading Kilosort config')
with open('../configs/kilosort_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

kilosort_path = config['base_path']
temp_data_path = config['temp_data_path']

# Select preprocessed ephys data
input_dir = askdirectory(title='Please select root of CatGT-processed directory', initialdir=config['data_path'])

epoch_name = os.listdir(input_dir)[0]
run_name = os.listdir(input_dir)[0][0:-3]
run_name = run_name.replace('catgt_', '')

# Get run info and number of probes/channels saved
dirnames = 1
n_probes = len(next(os.walk(os.path.join(input_dir, epoch_name)))[dirnames])
print('Recorded using {} probes'.format(n_probes))

for probe_id in range(n_probes)[2:4]:
    probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
    probe_path = os.path.join(input_dir, epoch_name, probe_folder)

    #binfile_name = [f for f in os.listdir(probe_path) if f.endswith('ap.bin')][0]
    #binfile_path = os.path.join(input_dir, epoch_name, probe_folder, binfile_name)

    #metafile_name = '{}_tcat.imec{}.ap.meta'.format(epoch_name.replace('catgt_', ''), probe_id)
    #metafile_path = os.path.join(input_dir, epoch_name, probe_folder, metafile_name)
    #meta_dict = readSGLX.readMeta(pathlib.Path(metafile_path))

    ## Add number of saved channels
    #print('-- IMEC probe', probe_id, ' #channels saved:', meta_dict['nSavedChans'])
    #n_saved_ch_probes[probe_id] = meta_dict['nSavedChans']

    # Start MATLAB engine
    eng = matlab.engine.start_matlab()
    eng.addpath(eng.genpath(r'C:\Users\bisi\Github\npy-matlab'), nargout=0)
    eng.cd(kilosort_path, nargout=0)

    # Run Kilosort for current probe
    print('- Running Kilosort for IMEC probe', probe_id)
    eng.run_main_kilosort(probe_path, temp_data_path, nargout=0)

    # Stop MATLAB engin
    eng.quit()
