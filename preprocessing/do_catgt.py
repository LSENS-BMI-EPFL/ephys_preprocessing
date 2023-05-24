#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: do_catgt.py
@time: 28/03/2022 15:07
"""

# Imports 
import os
import json
import pprint
import subprocess
import tkinter.filedialog as fdialog
import webbrowser
from pathlib import Path
import sys

from ephys_utils import flatten_list


# Modules
sys.path.append('C:\\Users\\bisi\\Github\\ecephys_spike_sorting\\ecephys_spike_sorting\\common')
from SGLXMetaToCoords import readMeta

# Load config
print('Loading CatGT config')
with open('../configs/catgt_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

# Select raw input ephys recording run
input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory', initialdir=config['data_path'])
session_name = [i for i in os.listdir( os.path.join(input_dir_mouse, 'Recording')) if 'AB' in i][0]
input_dir = os.path.join(input_dir_mouse, 'Recording', session_name, 'Ephys')

epoch_name = os.listdir(input_dir)[0]
epoch_number = epoch_name[-1]
run_name = os.listdir(input_dir)[0][0:-3]

print('Input data directory:', input_dir, run_name)


# Select output mouse directory
output_dir_mouse = fdialog.askdirectory(title='Please select mouse output directory', initialdir=config['save_path'])
output_dir = os.path.join(output_dir_mouse, 'Recording', session_name, 'Ephys')
Path(output_dir).mkdir(parents=True, exist_ok=True) # create output dir

print('Output directory', output_dir)

# Get run info and number of probes/channels saved
dirnames = 1 # takes first run i.e. gN_1st
n_probes = len(next(os.walk(os.path.join(input_dir, epoch_name)))[dirnames])
print('Recorded using {} probe(s)'.format(n_probes))
n_saved_ch_probes = dict.fromkeys(range(n_probes))

for probe_id in range(n_probes):
    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
    metafile_name = '{}_t0.imec{}.ap.meta'.format(epoch_name, probe_id)
    metafile_path = os.path.join(input_dir, epoch_name, probe_folder, metafile_name)

    meta_dict = readMeta(Path(metafile_path))

    # Add number of saved channels
    print('-- IMEC probe', probe_id, ' #channels saved:', meta_dict['nSavedChans'])
    n_saved_ch_probes[probe_id] = int(meta_dict['nSavedChans'])


# Write CatGT command line
command = ['CatGT',
           '-dir={}'.format(input_dir),
           '-run={}'.format(run_name),  #mouse name basically
           '-prb_fld', '-prb_miss_ok',  #assumes probe data saved in separate folders
           '-g={}'.format(epoch_number),                      #saved SGLX run not necessarily the first one (g-index)
           '-t=0,0','-t_miss_ok',      #assumes only one SGLX run
           #'-lf', '-ap',
           '-prb=0:5',                  #assumes at most 6 probes
           '-ni'
           ]

#for probe_id in range(n_probes): #old CatGT version
#    command.append(['-SY=2,{},{},6,500'.format(probe_id, n_saved_ch_probes[0]-1)])

command.append([
           #'-xd=2,{},{},6,500'.format(probe_id,  n_saved_ch_probes[0]-1),  #Commmented as performed by default
           '-xa=0,0,0,1,0,0',               #Square wave pulse from IMEC slot (alsone by default)
           '-xa=0,0,1,1,0,0',               #Trial start
           '-xa=0,0,2,0.5,1,0',             #Auditory stimulus (does not work)
           '-xa=0,0,3,0.5,1,0',             #Whisker stimulus
           '-xa=0,0,4,1,0,0',               #Valve opening
           #'-xa=0,0,5,1,0,0',               #Behaviour camera 0 frame times
           #'-xa=0,0,6,1,0,0',               #Behaviour camera 1 frame times
           #'-xa=0,0,7,0.005,0.010,0',       #Piezo lick sensor #TEST #second piezo?
           #'-gblcar',                        #global CAR
           '-dest={}'.format(output_dir),
           '-out_prb_fld'])             #saved in separate probe folders

print('CatGT command line will run:', list(flatten_list(command)))


# Run CatGT
print('Running CatGT...')
subprocess.run(list(flatten_list(command)), shell=True, cwd=config['base_path'])
print('CatGT done!')

# Open log file
print('Opening log file')
webbrowser.open(os.path.join(config['base_path'], 'CatGT.log'))

print('Finished CatGT.')
