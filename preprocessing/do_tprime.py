#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: do_tprime.py
@time: 01/04/2022 13:23
"""

# Imports
import os
import json
import pprint
import subprocess
import numpy as np
import tkinter.filedialog as fdialog  # gui
import webbrowser
from pathlib import Path
import sys
from collections.abc import Iterable

# Modules
sys.path.append('C:\\Users\\bisi\\Github\\ecephys_spike_sorting\\ecephys_spike_sorting\\common')
from SGLXMetaToCoords import readMeta  # read syncperiod

# Load config
print('Loading TPrime config')
with open('configs/tprime_config.json') as json_conf:
    config = json.load(json_conf, strict=False)
pprint.pprint(config)

# Select processed mouse data directory
output_dir_mouse = fdialog.askdirectory(title='Please select mouse output directory', initialdir=config['analysis_data_path'])
session_name = [i for i in os.listdir( os.path.join(output_dir_mouse, 'Recording')) if 'AB' in i][0]
output_dir = os.path.join(output_dir_mouse, 'Recording', session_name, 'Ephys')

catgt_epoch_name = [d for d in os.listdir(output_dir) if 'catgt' in d][0]
epoch_name = catgt_epoch_name[6:] #MOUSENAME_gX
run_name = catgt_epoch_name[6:-3] # MOUSENAME

# Get synchronization period
sglx_metafile_path = Path(output_dir, catgt_epoch_name, '{}_tcat.nidq.meta'.format(epoch_name))
sglx_meta_dict = readMeta(sglx_metafile_path)

if sglx_meta_dict['syncSourcePeriod'] is None:
    syncperiod = float(config['syncperiod'])
else:
    syncperiod = float(sglx_meta_dict['syncSourcePeriod'])

# Get number of probes
dirnames = 1
subfolder_list = next(os.walk(os.path.join(output_dir, catgt_epoch_name)))[dirnames]
probe_folder_list = [s for s in subfolder_list if 'imec' in s]
n_probes = len(probe_folder_list)
print('Recorded using {} probes'.format(n_probes))

# Probe recordings with spikes (to include)
valid_probes = []

## Convert spike times from samples to seconds
print('Convert spike times in seconds')
for probe_id in range(n_probes):
    print('-- IMEC probe {} spike times in seconds'.format(probe_id))

    # Get probe SR
    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
    metafile_name = '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id)
    apbin_metafile_path = os.path.join(output_dir, catgt_epoch_name, probe_folder, metafile_name)

    #Read probe meta information
    ap_meta_dict = readMeta(Path(apbin_metafile_path))
    imSampRate = float(ap_meta_dict['imSampRate'])  # probe-specific

    try:
        # Load spike times
        spike_times = np.load(os.path.join(output_dir, catgt_epoch_name, probe_folder, 'spike_times.npy'))
        # Convert into seconds
        np.save(os.path.join(output_dir, catgt_epoch_name, probe_folder, 'spike_times_sec.npy'), spike_times / imSampRate)
        valid_probes.append(probe_id)

    except FileNotFoundError as e:
        print('No spike times for IMEC probe {}: either spike sorting missing or bad recording'.format(probe_id))

## Write Tprime command line
print('Assembling TPrime command')
nidq_stream_idx = 10 #arbitrary index number

## Set path to reference alignment probe
# First check default reference probe is included
if config['default_tostream_probe'] not in valid_probes:
    default_tostream_probe = config['default_tostream_probe'] + 1 #else use next one
else:
    default_tostream_probe = config['default_tostream_probe']

path_ref_probe = os.path.join(output_dir, catgt_epoch_name, '{}_imec{}'.format(epoch_name,
                                                             default_tostream_probe))
ref_probe_edges_file = '{}_tcat.imec{}.ap.xd_{}_6_500.txt'.format(epoch_name,
                                                                  default_tostream_probe,
                                                                  int(ap_meta_dict['nSavedChans'])-1)

# Create output folder with aligned event times
path_dest = os.path.join(output_dir, catgt_epoch_name, 'sync_event_times')
Path(path_dest).mkdir(parents=True, exist_ok=True)


# Set reference streams
command = ['Tprime',
           '-syncperiod={}'.format(syncperiod),
           # arg: reference data stream edge times (IMEC 0)
           '-tostream={}'.format(os.path.join(path_ref_probe, ref_probe_edges_file)),
           # arg: stream index, sync pulse edge times
           '-fromstream={},{}'.format(nidq_stream_idx,
                                      os.path.join(output_dir, catgt_epoch_name, epoch_name+'_tcat.nidq.xa_0_0.txt'))
           ]

# Add edge times & spike times for included each probe
for probe_id in valid_probes:
    print('-- IMEC probe {} spike times sync arguments added'.format(probe_id))

    probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
    probe_folder_path = os.path.join(output_dir, catgt_epoch_name, probe_folder)

    probe_edgetime_files = [f for f in os.listdir(probe_folder_path) if 'ap.xd' in f]

    command.append('-fromstream={},{}'.format(probe_id,
                                              os.path.join(probe_folder_path,
                                                           probe_edgetime_files[0])))  # stream index, edge times (probe)
    command.append('-events={},{},{}'.format(probe_id,
                                             os.path.join(probe_folder_path, 'spike_times_sec.npy'),
                                             os.path.join(probe_folder_path, #save in original imec folder
                                                          '{}_imec{}_spike_times_sec_sync.npy'.format(epoch_name,probe_id))))  # stream index, spike times
    command.append('-events={},{},{}'.format(probe_id,
                                             os.path.join(probe_folder_path, 'spike_times_sec.npy'),
                                             os.path.join(path_dest, #save along other aligned event times
                                                          '{}_imec{}_spike_times_sec_sync.npy'.format(epoch_name,
                                                                                                           probe_id))))


# Add behaviour/video event times
output_path = Path(output_dir, catgt_epoch_name)
command.append([
    '-events={},{},{}'.format(nidq_stream_idx,
                              os.path.join(output_path, '{}_tcat.nidq.xa_1_0.txt'.format(epoch_name)),
                              os.path.join(path_dest, 'trial_start_times.txt')),
    '-events={},{},{}'.format(nidq_stream_idx,
                              os.path.join(output_path, '{}_tcat.nidq.xa_2_0.txt'.format(epoch_name)),
                              os.path.join(path_dest, 'auditory_stim_times.txt')),
    '-events={},{},{}'.format(nidq_stream_idx,
                              os.path.join(output_path, '{}_tcat.nidq.xa_3_0.txt'.format(epoch_name)),
                              os.path.join(path_dest, 'whisker_stim_times.txt')),
    '-events={},{},{}'.format(nidq_stream_idx,
                              os.path.join(output_path, '{}_tcat.nidq.xa_4_0.txt'.format(epoch_name)),
                              os.path.join(path_dest, 'valve_times.txt')),
    #'-events={},{},{}'.format(nidq_stream_idx, #TODO: replace by camrea strobe out copies, 1 and 2
    #                          os.path.join(output_path, '{}_tcat.nidq.xa_5_0.txt'.format(epoch_name)),
    #                          os.path.join(path_dest, 'camera_frame_times.txt')),
    #'-events={},{},{}'.format(nidq_stream_idx,
    #                          os.path.join(output_path, '{}_tcat.nidq.xa_6_0.txt'.format(epoch_name)),
    #                          os.path.join(path_dest, 'camera_arming_times.txt')),
    #'-events={},{},{}'.format(nidq_stream_idx,
    #                          os.path.join(output_path, '{}_tcat.nidq.xa_7_0.txt'.format(epoch_name)),
    #                          os.path.join(path_dest, 'piezo_licks.txt')) #TODO: remove, not working?

])

def flatten(l):
    """
    A function to flatten a list of list.
    :param l: A list containing lists.
    :return: Generator of the iterable.
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


print('Tprime command line will run:', list(flatten(command)))

# Run Tprime
print('Running TPrime...')
subprocess.run(list(flatten(command)), shell=True, cwd=config['base_path'])
print('TPrime done!')

# Open log file
print('Opening log file')
webbrowser.open(os.path.join(config['base_path'], 'Tprime.log'))

print('Finished Tprime.')

