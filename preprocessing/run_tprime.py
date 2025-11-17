#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_tprime.py
@time: 8/2/2023 8:02 PM
"""

# Imports
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import subprocess
import webbrowser
import pathlib
import numpy as np
from loguru import logger

from utils import readSGLX
from utils.ephys_utils import flatten_list


def main(input_dir, config):
    """
    Run TPrime on processed and spike-sorted/curated ephys data.
    This aligns task events and spike times to the same time base.
    :param input_dir:  path to CatGT processed ephys data
    :param config:  config dict
    :return:
    """

    catgt_epoch_name = pathlib.Path(input_dir).name

    # Create output folder with aligned event times
    path_dest = os.path.join(input_dir, 'sync_event_times')
    pathlib.Path(path_dest).mkdir(parents=True, exist_ok=True)

    # Get epoch name and run name
    epoch_name = catgt_epoch_name[6:]  # MOUSENAME_gX

    # Get synchronization period
    sglx_metafile_path = os.path.join(input_dir, '{}_tcat.nidq.meta'.format(epoch_name))
    sglx_meta_dict = readSGLX.readMeta(pathlib.Path(sglx_metafile_path))

    # Use specified syncperiod if available, otherwise use default
    if sglx_meta_dict['syncSourcePeriod'] is None:
        syncperiod = float(config['syncperiod'])
    else:
        syncperiod = float(sglx_meta_dict['syncSourcePeriod'])

    # Get number of probes
    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    n_probes = len(probe_folders)

    valid_probes = []
    for probe_id in range(n_probes):
        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        metafile_name = '{}_tcat_corrected.imec{}.ap.meta'.format(epoch_name, probe_id)
        apbin_metafile_path = os.path.join(input_dir, probe_folder, metafile_name)
        ap_meta_dict = readSGLX.readMeta(pathlib.Path(apbin_metafile_path))
        imSampRate = float(ap_meta_dict['imSampRate'])  # probe-specific

        try:
            logger.info('Converting IMEC probe {} spike times to seconds.'.format(probe_id))
            # Find any folders that contain kilosort in the name
            kilosort_folder = 'kilosort2' # TODO: if other KS version used, generalize this

            # Load spike times and convert in seconds
            spike_times = np.load(os.path.join(input_dir, probe_folder, kilosort_folder, 'spike_times.npy'))
            spike_times_sec = spike_times / imSampRate
            path_to_spikes = os.path.join(input_dir, probe_folder, kilosort_folder, 'spike_times_sec.npy')
            np.save(path_to_spikes, spike_times_sec)
            valid_probes.append(probe_id)

        # If no spike times, skip probe
        except FileNotFoundError as e:
            logger.warning('No spike times for IMEC probe {}: either spike sorting missing or invalid recording.'.format(probe_id))

    # Write TPrime command line
    nidq_stream_idx = 10  # arbitrary index number

    ## Set path to reference alignment probe
    # First check default reference probe is included in valid probes
    if config['default_tostream_probe'] not in valid_probes:
        default_tostream_probe = config['default_tostream_probe'] + 1  # else use next one
    else:
        default_tostream_probe = config['default_tostream_probe']

    path_ref_probe = os.path.join(input_dir, '{}_imec{}'.format(epoch_name, default_tostream_probe))
    ref_probe_edges_file = '{}_tcat.imec{}.ap.xd_{}_6_500.txt'.format(epoch_name,
                                                                      default_tostream_probe,
                                                                      int(ap_meta_dict['nSavedChans']) - 1)
    # Set reference streams
    command = ['Tprime',
               '-syncperiod={}'.format(syncperiod),                                         # arg: reference data stream edge times (IMEC 0)
               '-tostream={}'.format(os.path.join(path_ref_probe, ref_probe_edges_file)),   # arg: stream index, sync pulse edge times
               '-fromstream={},{}'.format(nidq_stream_idx,
                                          os.path.join(input_dir, epoch_name + '_tcat.nidq.xa_0_0.txt'))
               ]

    # Add edge times & spike times for included each probe
    for probe_id in valid_probes:
        logger.info('Adding IMEC probe {} spike times sync arguments'.format(probe_id))

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        probe_folder_path = os.path.join(input_dir, probe_folder)

        probe_edgetime_files = [f for f in os.listdir(probe_folder_path) if 'ap.xd' in f]

        command.append('-fromstream={},{}'.format(probe_id,
                                                  os.path.join(probe_folder_path,
                                                               probe_edgetime_files[
                                                                   0])))  # stream index, edge times (probe)

        command.append('-events={},{},{}'.format(probe_id,
                                                 os.path.join(probe_folder_path, kilosort_folder, 'spike_times_sec.npy'),
                                                 os.path.join(probe_folder_path,  # save in original imec folder
                                                              '{}_imec{}_spike_times_sec_sync.npy'.format(epoch_name,
                                                                                                          probe_id))))  # stream index, spike times
        command.append('-events={},{},{}'.format(probe_id,
                                                 os.path.join(probe_folder_path, kilosort_folder, 'spike_times_sec.npy'),
                                                 os.path.join(path_dest,  # save AGAIN along other aligned event times
                                                              '{}_imec{}_spike_times_sec_sync.npy'.format(epoch_name,
                                                                                                          probe_id))))

# Add behaviour and video frame times
    command.append([
        '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_1_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'trial_start_times.txt')),
        '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_2_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'auditory_stim_times.txt')),
        '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_3_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'whisker_stim_times.txt')),
        #'-events={},{},{}'.format(nidq_stream_idx,
        #                          os.path.join(input_dir, '{}_tcat.nidq.xa_4_0.txt'.format(epoch_name)),
        #                          os.path.join(path_dest, 'valve_times.txt')), #TODO: for AB mice
        #'-events={},{},{}'.format(nidq_stream_idx,
        #                          os.path.join(input_dir, '{}_tcat.nidq.xa_4_0.txt'.format(epoch_name)),
        #                          os.path.join(path_dest, 'context_transition_on.txt')), # TODO: for PB mice
        #'-events={},{},{}'.format(nidq_stream_idx,
        #                          os.path.join(input_dir, '{}_tcat.nidq.xia_4_0.txt'.format(epoch_name)),
        #                          os.path.join(path_dest, 'context_transition_off.txt')),  # TODO: for PB mice
         '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_5_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'cam0_frame_times.txt')),
         '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_6_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'cam1_frame_times.txt')),
         '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_7_0.txt'.format(epoch_name)), # Note: this works weirdly
                                  os.path.join(path_dest, 'piezo_licks.txt'))

    ])
    if epoch_name.startswith('AB'):
        command.append(['-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_4_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'valve_times.txt'))
                        ])
    elif epoch_name.startswith('PB'):
        command.append(['-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xa_4_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'context_transition_on.txt')),
                        '-events={},{},{}'.format(nidq_stream_idx,
                                  os.path.join(input_dir, '{}_tcat.nidq.xia_4_0.txt'.format(epoch_name)),
                                  os.path.join(path_dest, 'context_transition_off.txt'))
                        ])

    logger.info('TPrime command line will run: {}'.format(list(flatten_list(command))))

    logger.info('Running TPrime to align task events and spike times.')
    subprocess.run(list(flatten_list(command)), shell=True, cwd=config['tprime_path'])

    logger.info('Opening TPrime log file at: {}'.format(os.path.join(config['tprime_path'], 'Tprime.log')))
    webbrowser.open(os.path.join(config['tprime_path'], 'Tprime.log'))

    return

