#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: recording.py
@time: 10/3/2022 3:15 PM
"""


class Recording:
    # RECORDING Organizes behavioural and neural data.


    def __init__(self, mouse_id):
        """

        :param mouse_id: (int) mouse number
        """

        self.probe_list =