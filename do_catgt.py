#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: do_catgt.py
@time: 28/03/2022 15:07
"""

# Imports 
import os
import json
import subprocess
import tkinter.filedialog as fdialog #gui
import webbrowser
from pathlib import Path
import sys
from collections.abc import Iterable

# Modules
sys.path.append('C:\\Users\\bisi\\Github\\ecephys_spike_sorting\\ecephys_spike_sorting\\common')
from SGLXMetaToCoords import readMeta

# Set paths
with open('catgt_config.json') as json_conf :
    config = json.load(json_conf, strict=False)
for key,val in config.items():
    print(key,val)

# Select raw input ephys recording run
input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory', initialdir=config['data_path'])
input_dir = os.path.join(input_dir_mouse, 'Recording/Ephys')

epoch_name = os.listdir(input_dir)[0]
run_name = os.listdir(input_dir)[0][0:-3]

# Select output mouse directory
output_dir_mouse = fdialog.askdirectory(title='Please select mouse output directory', initialdir=config['save_path'])
output_dir = os.path.join(output_dir_mouse, 'Recording/Ephys')

# Get run info and number of channels saved
dirnames = 1
n_probes = len(next(os.walk(os.path.join(input_dir, epoch_name)))[dirnames]) #count number of probe folder
print('Recorded using {} probes'.format(n_probes))
n_saved_ch_probes = dict.fromkeys(range(n_probes))

for probe_id in range(n_probes):
    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
    metafile_name = '{}_t0.imec{}.ap.meta'.format(epoch_name, probe_id)
    metafile_path = os.path.join(input_dir, epoch_name, probe_folder, metafile_name)

    meta_dict = readMeta(Path(metafile_path))

    # Add number of saved channels
    print('-- IMEC probe', probe_id, ' #channels saved:', meta_dict['nSavedChans'])
    n_saved_ch_probes[probe_id] = meta_dict['nSavedChans']


# Write CatGT command line
command = ['CatGT',
           '-dir={}'.format(input_dir),
           '-run={}'.format(run_name),  #mouse name basically
           '-prb_fld', '-prb_miss_ok',  #assumes probe data saved in separate folders
           '-g=0',                      #assumes only one SGLX run
           '-t=0,0', '-t_miss_ok',      #assumes only one SGLX run
           '-lf', '-ap',
           '-prb=0:5',                  #assumes at most 6 probes
           '-ni'
           ]
for probe_id in range(n_probes):
    command.append(['-SY={},{},6,500'.format(probe_id,n_saved_ch_probes[0])])

command.append([
           '-XD=8,0,500',               #?edges?
           '-XA=1,1,0,0',               #FILL IN
           '-XA=2,5,0,0',               #FILL IN
           '-XA=3,5,0,0',               #FILL IN
           '-XA=4,5,0,0',               #FILL IN
           '-XA=5,5,0,0',               #FILL IN
           '-XA=6,5,0,0',               #FILL IN
           '-XA=7,5,0,0',               #FILL IN
           '-gblcar',                   #global CAR
           '-dest={}'.format(output_dir),
           '-out_prb_fld'])             #saved in separate probe folders

def flatten(l):
    """
    A function to flatten a list of list.
    :param l: A list containing lists.
    :return:
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


print('CatGT command line will run:', list(flatten(command)))

# Run CatGT
subprocess.run(list(flatten(command)), shell=True, cwd=config['base_path'])

# Open log file
print('Opening log file')
webbrowser.open(os.path.join(config['base_path'], 'CatGT.log'))

print('Finished CatGT.')