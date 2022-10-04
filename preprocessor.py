#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: preprocessor.py
@time: 10/3/2022 10:09 AM
"""

# TODO : implement class
# Get path to recording folder in Npx_Data
# Read SGLX information using readSGLX
# Do CatGT
# Do KS
# Do TPrime
# Do C_Waves
# Format data, create files
# Implement method, eg coil artefact correction
# More ?

# Imports
import os, sys
import json
import glob
import numpy as np
import pandas as pd
import pprint as pp

class Preprocessor:
    # PREPROCESSOR Preprocessing of ephys data after recordings.

    def __init__(self):