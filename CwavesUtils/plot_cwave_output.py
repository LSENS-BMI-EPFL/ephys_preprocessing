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
    - mean waveforms,
    - waveform in probe geometry
    - C_waves waveform SNR information.

    :param m_name: Mouse name.
    :return:
    """

    # Load cluster SNR C_waves output
    cluster_snr = np.load(os.path.join(PATH_ANALYSIS,
                                       '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/cwaves/cluster_snr.npy'.format(m_name)))

    # Plot cluster SNRs distribution, spikes per pk. chan. distribution
    fig = plot_cluster_snr_hist(cluster_snr)
    fig.savefig(fname=os.path.join(PATH_ANALYSIS,
                                   '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/cwaves/cluster_snr_hist.png'.format(m_name)),
                bbox_inches='tight')

    fig = plot_cluster_spks_per_pkch_hist(cluster_snr)
    fig.savefig(fname=os.path.join(PATH_ANALYSIS,
                                   '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/cwaves/cluster_spks_pk_ch_hist.png'.format(m_name)),
                bbox_inches='tight')



    # Load mean waveforms C_waves output
    mean_wfs = np.load(os.path.join(PATH_ANALYSIS,
                                    '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/cwaves/mean_waveforms.npy'.format(m_name)))

    # Load cluster info phy output
    clus_info = pd.read_csv(os.path.join(PATH_ANALYSIS,
                                         '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/ks25/cluster_info.tsv'.format(m_name)),
                                         sep='\t')
    # Set index s.t. it matches cluster_id range (necessary for C_waves)
    clus_info.set_index(keys='cluster_id', drop=False, inplace=True)
    clus_info = clus_info.reindex(range(np.max(clus_info.cluster_id) + 1), fill_value=0, copy=True)

    # Create figure directory
    folder_path = os.path.join(PATH_ANALYSIS,
                               '{0}/Recording/Ephys/{0}_g0/{0}_g0_imec0/cwaves/pk_ch_mean_wf/'.format(m_name))
    Path(folder_path).mkdir(parents=True, exist_ok=True)

    # Plot all first 100 clusters
    for cluster_id in np.unique(clus_info.cluster_id)[:100]:

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
