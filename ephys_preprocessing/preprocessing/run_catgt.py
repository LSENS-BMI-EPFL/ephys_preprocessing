#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_catgt.py
@time: 8/2/2023 4:39 PM
"""


# Imports
import os
import sys
import subprocess
from ephys_preprocessing.utils.ephys_utils import flatten_list
import webbrowser
from loguru import logger

def main(input_dir, output_dir, config):
    """
    Run CatGT on raw ephys data and save to output directory.
    :param input_dir: path to raw ephys data
    :param config: config dict
    :return:
    """

    # Get epoch number and run name
    epoch_name = os.listdir(input_dir)[0]
    epoch_number = epoch_name[-1]
    run_name = os.listdir(input_dir)[0][0:-3]

    # Write CatGT command line
    if sys.platform.startswith('win'):
        catGTexe_fullpath = 'CatGT'
    elif sys.platform.startswith('linux'):
        catGTexe_fullpath = config['catgt_path']
        catGTexe_fullpath = catGTexe_fullpath.replace('\\', '/') + "/runit.sh"

    command = [catGTexe_fullpath,
               '-dir={}'.format(input_dir),
               '-run={}'.format(run_name),
               '-prb_fld',
               '-prb_miss_ok',
               '-g={}'.format(epoch_number),
               '-t=0,0',
               '-t_miss_ok',
               '-startsecs=0.0',
               '-maxsecs=5580', # remove and write down for mouse in SLIMS
               '-ni',
               '-lf',
               '-ap',
               '-prb=0:5',
               '-xa=0,0,0,1,0,0',               # Square wave pulse from IMEC slot (on by default)
               '-xa=0,0,1,4,0,0',               # Trial start
               '-xa=0,0,2,1,1,0',               # Auditory stimulus (does not work)
               '-xa=0,0,3,1,1,0',               # Whisker stimulus
               '-xa=0,0,4,2,0,0',               # Valve opening
               '-xa=0,0,5,2,0,0',               # Behaviour camera 0 frame times
               '-xa=0,0,6,2,0,0',               # Behaviour camera 1 frame times
               '-xa=0,0,7,0.005,0.010,0',       # Piezo lick sensor
               '-gblcar',                       # global common median referencing (default), never applied to LFP
               '-dest={}'.format(output_dir),
               '-out_prb_fld'
               ]

    logger.info('CatGT command line will run: {}'.format(list(flatten_list(command))))

    logger.info('Running CatGT on {}.'.format(epoch_name))
    if sys.platform.startswith('win'):
        subprocess.run(list(flatten_list(command)), shell=True, cwd=config['catgt_path'])
    elif sys.platform.startswith('linux'):
        subprocess.run(list(flatten_list(command)), cwd=config['catgt_path'])

    logger.info('Opening CatGT log file at: {}'.format(os.path.join(config['catgt_path'], 'CatGT.log')))
    webbrowser.open(os.path.join(config['catgt_path'], 'CatGT.log'))

    return
