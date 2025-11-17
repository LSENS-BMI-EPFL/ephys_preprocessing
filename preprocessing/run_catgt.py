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
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import subprocess
from utils.ephys_utils import flatten_list
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
    epoch_name = [f for f in os.listdir(input_dir) if '_g' in f][0]
    epoch_number = epoch_name[-1]
    run_name = epoch_name[0:-3]

    # Write CatGT command line
    command = ['CatGT',
               '-dir={}'.format(input_dir),
               '-run={}'.format(run_name),
               '-prb_fld',
               '-prb_miss_ok',
               '-g={}'.format(epoch_number),
               '-t=0,0',
               '-t_miss_ok',
               '-startsecs=0.0',
               #'-maxsecs=2998.0',    # TODO: remove (for supercat)
               #'-pass1_force_ni_ob_bin',# TODO: remove (for supercat)
               '-ni',
               '-lf',
               '-ap',
               '-prb=0:5',
               '-xa=0,0,0,1,0,0',               # Square wave pulse from IMEC slot (on by default)
               '-xa=0,0,1,4,0,0',               # Trial start
               '-xa=0,0,2,1,1,0',               # Auditory stimulus (does not work)
               '-xa=0,0,3,1,1,0',               # Whisker stimulus
               #'-xa=0,0,4,2,0,0',               # Valve opening #TODO: for AB mice
               #'-xa=0,0,4,2,0,0',               # Context transition TTL epoch start #TODO: PB mice
               #'-xia=0,0,4,2,0,0',               # Context transition TTL epoch end#TODO: PB mice
               '-xa=0,0,5,2,0,0',               # Behaviour camera 0 frame times
               '-xa=0,0,6,2,0,0',               # Behaviour camera 1 frame times
               '-xa=0,0,7,0.005,0.010,0',       # Piezo lick sensor
               '-gblcar',                       # global common median referencing (default), never applied to LFP
               '-dest={}'.format(output_dir),
               '-out_prb_fld'
               ]
    if epoch_name.startswith('AB'):
        command.append(['-xa=0,0,4,2,0,0']) # valve opening times
    elif epoch_name.startswith('PB'):
        command.append(['-xa=0,0,4,2,0,0',  # context transition TTL epoch start
                        '-xia=0,0,4,2,0,0']) # context transition TTL epoch end



    logger.info('CatGT command line will run: {}'.format(list(flatten_list(command))))

    logger.info('Running CatGT on {}.'.format(epoch_name))
    subprocess.run(list(flatten_list(command)), shell=True, cwd=config['catgt_path'])

    logger.info('Opening CatGT log file at: {}'.format(os.path.join(config['catgt_path'], 'CatGT.log')))
    webbrowser.open(os.path.join(config['catgt_path'], 'CatGT.log'))

    return
