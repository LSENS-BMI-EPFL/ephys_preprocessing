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
import subprocess
import webbrowser
import pathlib
import numpy as np
from loguru import logger

from ephys_preprocessing.utils import readSGLX
from ephys_preprocessing.utils.ephys_utils import flatten_list


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
        metafile_name = '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id)
        apbin_metafile_path = os.path.join(input_dir, probe_folder, metafile_name)
        ap_meta_dict = readSGLX.readMeta(pathlib.Path(apbin_metafile_path))
        imSampRate = float(ap_meta_dict['imSampRate'])  # probe-specific

        kilosort_folders = (pathlib.Path(input_dir) / probe_folder).glob('kilosort*')
        for kilosort_folder in kilosort_folders:

            try:
                logger.info('Converting IMEC probe {} spike times to seconds for {}.'.format(probe_id, kilosort_folder.name))

                # Load spike times and convert in seconds
                spike_times = np.load(kilosort_folder / 'sorter_output' / 'spike_times.npy')
                spike_times_sec = spike_times / imSampRate
                path_to_spikes = kilosort_folder / 'sorter_output' / 'spike_times_sec.npy'
                np.save(path_to_spikes, spike_times_sec)
                valid_probes.append(probe_id)

            # If no spike times, skip probe
            except FileNotFoundError as e:
                logger.warning('No spike times for IMEC probe {} for {}: either spike sorting missing or invalid recording.'.format(probe_id, kilosort_folder.name))

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
    if sys.platform.startswith('win'):
        Tprime_fullpath = 'Tprime'
        shell = True
    elif sys.platform.startswith('linux'):
        Tprime_fullpath = config['tprime_path']
        Tprime_fullpath = Tprime_fullpath.replace('\\', '/') + "/runit.sh"
        shell = False
    else:
        raise NotImplementedError('OS not recognised')

    command = [Tprime_fullpath,
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

        kilosort_folders = [ks.name for ks in pathlib.Path(probe_folder_path).glob('kilosort*')]
        for kilosort_folder in kilosort_folders:
            command.append('-events={},{},{}'.format(probe_id,
                                                 os.path.join(probe_folder_path, kilosort_folder, 'sorter_output', 'spike_times_sec.npy'),
                                                 os.path.join(probe_folder_path,  # save in original imec folder
                                                              '{}_imec{}_{}_spike_times_sec_sync.npy'.format(epoch_name,
                                                                                                          probe_id,
                                                                                                          kilosort_folder))))  # stream index, spike times
            command.append('-events={},{},{}'.format(probe_id,
                                                 os.path.join(probe_folder_path, kilosort_folder, 'sorter_output', 'spike_times_sec.npy'),
                                                 os.path.join(path_dest,  # save AGAIN along other aligned event times
                                                              '{}_imec{}_{}_spike_times_sec_sync.npy'.format(epoch_name,
                                                                                                          probe_id,
                                                                                                          kilosort_folder))))

    # Define channel mappings for different setups
    SETUP_CONFIGS = {
        543: {  # Myri's setup
            'name': 'Myri setup',
            'channels': {
                'trial_start_times': 1,
                'whisker_stim_times': 2,
                'piezo_licks': 3,
                # 'valve_times': 4,
                'context_transition_on': 4,
                'context_transition_off': 4,  # Using xia_4_0
                'cam0_frame_times': 5,
                'cam1_frame_times': 6,
                'auditory_stim_times': 7
            }
        },
        'default': {  # Axel's setup
            'name': 'Axel setup',
            'channels': {
                'trial_start_times': 1,
                'auditory_stim_times': 2,
                'whisker_stim_times': 3,
                # 'valve_times': 4,
                'context_transition_on': 4,
                'context_transition_off': 4,  # Using xia_4_0
                'cam0_frame_times': 5,
                'cam1_frame_times': 6,
                'piezo_licks': 7
            }
        }
    }

    # Get setup configuration
    setup_sn = float(ap_meta_dict['imDatBsc_sn'])
    setup_config = SETUP_CONFIGS.get(setup_sn, SETUP_CONFIGS['default'])
    print(f"{setup_config['name']}")

    # Build event commands
    event_commands = []
    for event_name, channel in setup_config['channels'].items():
        # Special handling for context transition off which uses xia instead of xa
        if event_name == 'context_transition_off':
            channel_file = f'{epoch_name}_tcat.nidq.xia_{channel}_0.txt'
        else:
            channel_file = f'{epoch_name}_tcat.nidq.xa_{channel}_0.txt'
            
        event_commands.append(
            '-events={},{},{}'.format(
                nidq_stream_idx,
                os.path.join(input_dir, channel_file),
                os.path.join(path_dest, f'{event_name}.txt')
            )
        )

    command.extend(event_commands)

    logger.info('TPrime command line will run: {}'.format(list(flatten_list(command))))

    logger.info('Running TPrime to align task events and spike times.')
    subprocess.run(list(flatten_list(command)), shell=shell, cwd=config['tprime_path'])

    # logger.info('Opening TPrime log file at: {}'.format(os.path.join(config['tprime_path'], 'Tprime.log')))
    # webbrowser.open(os.path.join(config['tprime_path'], 'Tprime.log'))

    return

