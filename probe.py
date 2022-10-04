#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: probe.py
@time: 10/3/2022 2:50 PM
"""
# Imports
import sys, glob

# Modules
sys.path.append(r'C:\Users\bisi\Github\SpikeGLX_Datafile_Tools\SpikeGLX_Datafile_Tools\Python\DemoReadSGLXData')
from DemoReadSGLXData import readSGLX

class Probe:
    # PROBE Probe class for recorded data and experiments.

    def __init__(self, mouse_id, imec_id):
        """
        :param mouse_id: (int) mouse number
        :param imec_id: (int) imec probe identifier (starts at 0)
        """
        # General attributes
        self.rec_mouse = 'AB' + str(0) + str(mouse_id)
        self.rec_folder_path = None

        self.ap_bin_path = None
        self.ap_meta_path = None
        self.lf_bin_path = None
        self.lf_meta_path = None

        # Neural data
        self.cluster_info = []
        self.spikes_times = []
        self.spike_clusters = []

        # Anatomy
        self.zaxis_depth = None
        self.atlas_coords = {}
        self.channel_areas = []


    def get_sglx_information(self):
        """
        """
