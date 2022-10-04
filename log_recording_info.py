#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: log_recording_info.py
@time: 30/05/2022 18:28
"""

# Imports
import os
import json
from pathlib import Path
from datetime import datetime

# Paths
SAVED_DATA_PATH = r'D:\Npx_Data'

# Get information
mouse_name = input("Enter mouse name")
rec_date = datetime.today().strftime('%d%m%Y')

# Init. dict
rec_info = dict()

# Get general information
rec_info['mouse_name'] = [mouse_name]
value = input("Did you record today (date)?{}".format(p))
if value == 'yes':
    rec_info['date'] = [rec_date]
elif value == 'no':
    date = input("Enter recording date (dd.mm.yyyy):".format(p))
    rec_info['date'] = [rec_date]


# Get number of probes used
value = int(input("Enter number of probes used:"))
rec_info['n_probes'] = [value]

# Fill in info. for each probe
for p in range(value):
    rec_info_probe = dict()
    value = input("Enter cortical surface area name for IMEC probe {} (abbreviation):".format(p))
    rec_info_probe['target_area'] = [value]

    value = input("Enter insertion depth for IMEC probe {} (microns):".format(p))
    rec_info_probe['depth'] = [value]

    value = input("Enter targeting method for IMEC probe {} (ios/coord):".format(p))
    rec_info_probe['target_method'] = [value]

    value = input("Enter elevation angle (rel. to horizontal) for IMEC probe {} (read mount angle):".format(p))
    rec_info_probe['elevation_angle'] = [value]

    value = input("Enter azimuth angle (rel. to nose-tail line) for IMEC probe {} (read circular base: right positive, left negative):".format(p))
    rec_info_probe['azimuth_angle'] = [value]

    rec_info['imec{}'.format(p)] = rec_info_probe

base = Path(SAVED_DATA_PATH)
mouse_run_folders = [f for f in os.listdir(base) if mouse_name in f]
print('SpikeGLX runs recorded:', mouse_run_folders)

# Save recorded info as .log file in each run folder
for f in mouse_run_folders:
    jsonfile_path = base / f / '{}_{}_rec_info_dict.json'.format(mouse_name, rec_date)
    jsonfile_path.write_text(json.dumps(rec_info))
    print('Saved recording info as .json file in ', jsonfile_path)