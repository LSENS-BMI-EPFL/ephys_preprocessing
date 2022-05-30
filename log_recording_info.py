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

mouse_name = input("Enter mouse name")
rec_date = datetime.today().strftime('%d%m%Y')

# Init. dict
rec_info = dict()

rec_info['mouse_name'] = [mouse_name]
rec_info['date'] = [rec_date]

# Get number of probes used
value = int(input("Enter number of probes used:"))

rec_info['n_areas'] = [value]

# Fill in info. for each probe
for p in range(value):
    rec_info_probe = dict()
    value = input("Enter surface area for IMEC probe {}:".format(p))
    rec_info_probe['target_area'] = [value]

    value = input("Enter insertion depth (microns) for IMEC probe {}:".format(p))
    rec_info_probe['depth'] = [value]

    value = input("Enter localization method for IMEC probe {} (ios/coord):".format(p))
    rec_info_probe['loc_method'] = [value]

    value = input("Enter vertical angle method for IMEC probe {}: (mount)".format(p))
    rec_info_probe['vert_angle'] = [value]

    value = input("Enter mediolateral angle method for IMEC probe {}: (circular base)".format(p))
    rec_info_probe['ml_angle'] = [value]


    rec_info['imec{}'.format(p)] = rec_info_probe

# Save as json
base = Path(r'D:\Npx_Data')
jsonfile_path = base / '{}_{}_rec_info_dict.log'.format(mouse_name, rec_date)
jsonfile_path.write_text(json.dumps(rec_info))
