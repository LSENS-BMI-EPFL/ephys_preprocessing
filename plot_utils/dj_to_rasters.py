#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: dj_to_rasters.py
@time: 10.05.2023 10:42
"""

# Imports
import os
import datajoint as dj
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from pathlib import Path

from utils import remove_top_right_frame
from lsens_datajoint import lsens_base
from ephys_datajoint import ephys


def make_raster_plot(cluster_key, trial_start_times, trial_ids_dict, trial_type_dict):

    # Get cluster data from datajoint
    cluster_data_df = pd.DataFrame(cluster_key)

    # Get spikes and cluster information
    c_spk_times = np.asarray(cluster_data_df['spike_times'].values)
    cluster_id = cluster_data_df['cluster_id'].values[0]
    depth = cluster_data_df['depth'].values[0]  # distance relative to probe tip
    # area = cluster_data_df['area'].values[0] #TODO :
    title_str = 'id {}, {:.0f}$\mu m$'.format(cluster_id, depth)

    # Plot settings

    pre_event_win = 0.02  # in sec
    post_event_win = 0.05  # in sec
    time_ticks = np.arange(-pre_event_win, post_event_win, 0.01)
    time_ticks_labels = np.arange(-pre_event_win, post_event_win, 0.01) * 1e3  # in msec
    time_ticks_labels = np.round(time_ticks_labels).astype(int)  # format as int
    n_trials = len(trial_start_times)
    trial_ticks = np.arange(0, n_trials, 100)
    ft_size = 13

    #line_prop = dict(joinstyle='miter')

    # Make figure
    fig, ax = plt.subplots(1, 1, figsize=(2, 3), dpi=200)
    remove_top_right_frame(ax)
    ax.set_ylabel('Trials', fontsize=ft_size - 4)
    ax.set_xlabel('Time [ms]', fontsize=ft_size - 4)
    ax.set_yticks(ticks=trial_ticks, labels=trial_ticks, fontsize=ft_size - 5)
    ax.set_xticks(ticks=time_ticks, labels=time_ticks_labels, fontsize=ft_size - 5)
    ax.set_title(title_str, fontsize=ft_size / 2)
    ax.tick_params(axis='both', which='major', labelsize=ft_size / 2)

    trial_type_delimiters = []
    # Iterate over trial types
    for idx, (t_type, t_color) in enumerate(trial_type_dict.items()):
        trial_type_starts = trial_start_times[trial_ids_dict[t_type]]
        trial_type_delimiters.append(len(trial_type_starts))

        c_spks_aligned = []
        # Iterate over ordered trials
        for t_time in trial_type_starts:
            start_sec = t_time - pre_event_win
            end_sec = t_time + post_event_win
            c_spk_times_win = c_spk_times[(np.where((c_spk_times >= start_sec) & (c_spk_times <= end_sec)))] - t_time
            c_spks_aligned.append(c_spk_times_win)

        # Add raster per trial type
        if idx == 0:
            offset_start = 0
        else:
            offset_start = np.cumsum(trial_type_delimiters)[idx - 1]
        offset_end = np.cumsum(trial_type_delimiters)[idx]
        ax.eventplot(positions=c_spks_aligned,
                     lineoffsets=np.arange(offset_start, offset_end),
                     linewidths=1,
                     linelengths=6,
                     colors=[t_color] * len(c_spks_aligned),
                     linestyles='dotted',
                     )

    ax.axvline(x=0, lw=0.5, ls='--', c='k', zorder=0)
    ax.set_ylim(0, n_trials)
    ax.set_xlim(-pre_event_win, post_event_win)

    return fig


def dj_to_rasters(m_name):
    # Get mouse session data
    mouse_key = 'mouse_name="{}"'.format(m_name)
    sess_data = (ephys.EphysSession() * lsens_base.BehaviourTrial.Outcome() & mouse_key).fetch()
    sess_data_df = pd.DataFrame(sess_data)

    path_processed = (ephys.EphysSession() & mouse_key).fetch('path_to_rec_folder')[0]
    path_rasters = os.path.join(os.path.dirname(path_processed), 'rasters')
    Path(path_rasters).mkdir(parents=True, exist_ok=True)

    # Get trial types
    trial_types = ['cr', 'fa', 'wh', 'wm', 'ah', 'am']
    trial_colors = ['dimgray', 'k', 'forestgreen', 'crimson', 'mediumblue', 'lightblue']
    trial_type_dict_default = OrderedDict(list(zip(trial_types, trial_colors)))
    existing_trial_types = [c for c in trial_types if (1 in sess_data_df[c].values)]
    trial_type_dict = OrderedDict(list(zip(existing_trial_types,
                                           [trial_type_dict_default[t_type] for t_type in existing_trial_types])))

    # Get trial starts, ordred
    trial_start_times = \
        (ephys.EphysSessionTimestamps() & mouse_key & 'ts_name="trial_start_times"').fetch('ts_array')[0]

    trial_ids_dict = {}
    for t_type in trial_type_dict.keys():
        trial_ids_dict[t_type] = sess_data_df[sess_data_df[t_type] == 1].index.values

    for cluster_key in (ephys.EphysCluster.SpikeTimes() * ephys.EphysCluster() & mouse_key):

        # Make plot
        fig = make_raster_plot(cluster_key, trial_start_times, trial_ids_dict, trial_type_dict)

        # Save figure
        probe_id = cluster_key['probe_id']
        cluster_id = cluster_key['cluster_id']
        fname = '{}_probe{}_unit{}.png'.format(m_name.lower(), probe_id, cluster_id)
        fig.savefig(os.path.join(path_rasters, fname), format='png', dpi='figure', bbox_inches='tight')

    return


if __name__ == '__main__':
    # Connect to datajoint
    dj.config["enable_python_native_blobs"] = True
    dj.config['database.host'] = 'datajoint.epfl.ch'  # 'localhost:3306'
    dj.config['database.user'] = 'bisi'  # enter your username
    dj.config['database.password'] = 'C4XFho2E26rlFMR'  # enter your password
    ConnMessage = dj.conn()

    # Get args
    parser = argparse.ArgumentParser()
    parser.add_argument('--m_name', type=str, nargs='?', default='ABXXX', required=False)
    args = parser.parse_args()

    # For every neuron
    dj_to_rasters(args.m_name)
