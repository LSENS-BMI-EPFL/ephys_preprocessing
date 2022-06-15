#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: plot_cwave_output.py
@time: 17/03/2022 14:51
"""

# Imports
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter.filedialog as fdialog
from pathlib import Path

# Paths
PATH_ANALYSIS = 'M:/analysis/Axel_Bisi/mice_data/'

def plot_cluster_snr_hist(cluster_snr_cwave):
    """
    Plot histogram SNR for all clusters from C_waves output.

    :param cluster_snr_cwave: Cluster SNR array C_Waves output.
    :return:
    """

    fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=200)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.hist(cluster_snr_cwave[:, 0], bins=20, ec='w')
    ax.set_xlim(0) #negative SNRs are set for excluded clusters (during phy stage)
    ax.set_xlabel('SNR')
    ax.set_ylabel('Cluster count')

    return fig

def plot_cluster_spks_per_pkch_hist(cluster_snr_cwave):
    """
    Plot histogram of number of spikes per peak channel for all clusters from C_waves output.

    :param cluster_snr_cwave: Cluster SNR array C_Waves output.
    :return:
    """

    fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=200)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.hist(cluster_snr_cwave[:, 1], bins=20, ec='w')
    ax.set_xlabel('Number of spikes in peak channel')
    ax.set_ylabel('Cluster count')

    return fig

def plot_mean_waveform_pk_ch(mean_wf_cwave, clus_info_df, cluster_id):
    """
    Plot mean waveform of cluster at peak channel.

    :param mean_wf_cwave: Mean waveform C_waves output.
    :param clus_info_df: Cluster info dataframe KS output.
    :param cluster_id: Cluster index.
    :return:
    """

    pk_chan_id = clus_info_df.loc[cluster_id, 'ch']  # get cluster peak channel
    cluster_label = clus_info_df.loc[cluster_id, 'KSLabel']  # KSlabel before phy

    fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=100)
    ax.set_title('cluster: {}, KSLabel: {}'.format(cluster_id, cluster_label))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sns.lineplot(ax=ax, x=range(82), y=mean_wf_cwave[cluster_id, pk_chan_id, :], lw=3)
    ax.set_ylabel(r'$\mu$V')
    ax.set_xlabel('Time [samples]') #30kHz SR

    return fig

def plot_mean_waveform_probe(mean_wf_cwave, clus_info_df, cluster_id):
    """
    Plot mean waveform of cluster at few channels around peak channel.

    :param mean_wf_cwave: Mean waveform C_waves output.
    :param clus_info_df: Cluster info dataframe KS output.
    :param cluster_id: Cluster index.
    :return:
    """

    pk_chan_id = clus_info_df.loc[cluster_id, 'ch']  # get cluster peak channel
    cluster_label = clus_info_df.loc[cluster_id, 'KSLabel'] # KSlabel before phy

    n_chan_show = 16  # number of probe channels to show

    fig, axs = plt.subplots(int(n_chan_show / 2), 2, figsize=(5, 8), dpi=100, sharex=True, sharey=True)
    plt.suptitle('pk chan: {}, KSlabel: {}'.format(pk_chan_id, cluster_label))

    # Array probe channels s.t. even ch.-> right column & odd ch.-> left column
    if pk_chan_id%2==0:
        channels_arr_ids = np.arange(pk_chan_id - int(n_chan_show / 2), pk_chan_id + int(n_chan_show / 2), 1)
    else:
        channels_arr_ids = np.arange(pk_chan_id+1 - int(n_chan_show / 2), pk_chan_id+1 + int(n_chan_show / 2), 1)
    for idx, (ch_idx, ax) in enumerate(zip(channels_arr_ids, axs.flat)):

        ax.set_title('ch {}'.format(ch_idx))
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().set_ticks([])

        if ch_idx % 2 == 1:
            ax.spines['left'].set_visible(False)
            ax.yaxis.set_visible(False)
        if idx % 4 == 0 or idx % 4 == 1:
            ax_pos = ax.get_position()
            ax_pos.x0 = ax_pos.x0 + 0.05
            ax.set_position(ax_pos, which='original')
        try:
            sns.lineplot(ax=ax,
                        x=range(82),
                        y=mean_wf_cwave[cluster_id, ch_idx, :],
                        lw=3)
        except IndexError as err:
            print(err, ' peak channel close to channel max range -- ', cluster_id)

    return fig


def plot_cwave_output(m_name):
    """
    Plot output from C_waves:
    - mean cluster waveforms,
    - waveform in probe geometry layout,
    - C_waves waveform SNR information.

    :param m_name: Mouse name.
    :return:
    """

    #Select mouse directory to generate plots from
    input_dir_mouse = fdialog.askdirectory(title='Please select raw recording directory',
                                           initialdir=PATH_ANALYSIS)
    mouse_dir = os.path.join(input_dir_mouse, 'Recording/Ephys')
    epoch_name = os.listdir(mouse_dir)[0]  # for raw data
    catgt_epoch_name = 'catgt_{}'.format(epoch_name)  # for CatGT processed data

    #Count number of probe recordings
    dirnames=1
    if Path(mouse_dir, catgt_epoch_name).is_dir():
        path_cwave_output = Path(mouse_dir, catgt_epoch_name)
    elif Path(mouse_dir, epoch_name).is_dir():
        path_cwave_output = Path(mouse_dir, epoch_name)
    else:
        print('Neural recordings not available in', PATH_ANALYSIS)

    n_probes = len(next(os.walk(path_cwave_output))[dirnames])

    #Plot for each probe
    for probe_id in range(n_probes):
        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        path_cwave_probe = os.path.join(path_cwave_output, probe_folder, 'cwaves')
        path_kilosort_probe = os.path.join(path_cwave_output, probe_folder, 'ks25')
        print('Plotting for IMEC probe at {}'.format(os.path.join(path_cwave_output, probe_folder)))

        ## Load cluster SNR C_waves output
        cluster_snr = np.load(os.path.join(path_cwave_probe, 'cluster_snr.npy'))

        # Plot cluster SNRs distribution
        fig = plot_cluster_snr_hist(cluster_snr)
        fig.savefig(fname=os.path.join(path_cwave_probe, 'cluster_snr_hist.png'), bbox_inches='tight')

        # Plot distribution of spikes per pk. chan
        fig = plot_cluster_spks_per_pkch_hist(cluster_snr)
        fig.savefig(fname=os.path.join(path_cwave_probe, 'cluster_spks_pk_ch_hist.png'), bbox_inches='tight')


        ## Load mean waveforms C_waves output
        mean_wfs = np.load(os.path.join(path_cwave_probe, 'mean_waveforms.npy'))

        # Load cluster info phy output
        clus_info = pd.read_csv(os.path.join(path_kilosort_probe, 'cluster_info.tsv'), sep='\t') #<- assumes Phy performed

        # Set index s.t. it matches cluster_id range (necessary for C_waves)
        clus_info.set_index(keys='cluster_id', drop=False, inplace=True)
        clus_info = clus_info.reindex(range(np.max(clus_info.cluster_id) + 1), fill_value=0, copy=True)

        # Create figure directory
        folder_path = os.path.join(path_cwave_probe, 'pk_ch_mean_wf')
        Path(folder_path).mkdir(parents=True, exist_ok=True)

        # Plot all clusters
        for cluster_id in np.unique(clus_info.cluster_id):

            # Plot mean waveform at peak channel
            fig = plot_mean_waveform_pk_ch(mean_wf_cwave=mean_wfs, clus_info_df=clus_info, cluster_id=cluster_id)
            fig.savefig(fname=os.path.join(folder_path, 'cluster{}_pk_ch_mean_wf.png'.format(cluster_id)),
                        bbox_inches='tight')

            # Plot mean waveform per cluster in probe
            fig = plot_mean_waveform_probe(mean_wf_cwave=mean_wfs, clus_info_df=clus_info, cluster_id=cluster_id)
            fig.savefig(fname=os.path.join(folder_path, 'cluster{}_probe_mean_wf.png'.format(cluster_id)),
                        bbox_inches='tight')

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--m_name', type=str, nargs='?', default='ABXXX', required=True)
    args = parser.parse_args()

    plot_cwave_output(args.m_name)
