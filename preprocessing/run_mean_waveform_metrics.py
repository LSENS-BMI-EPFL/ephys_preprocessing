#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_mean_waveform_metrics.py
@time: 10/24/2023 4:17 PM
"""

# Imports
import os
from utils import readSGLX
import numpy as np
import pathlib
import pandas as pd
from utils.waveform_metrics_utils import calculate_waveform_metrics_from_avg


def main(input_dir):
    """
    Run mean waveform metrics on preprocessed spike data.
    This computes metrics for each cluster waveforms after C_Waves.
    :param input_dir: path to CatGT preprocessed data
    :param config: config dict
    :return:
    """

    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Perform computations for each probe separately
    for probe_id in probe_ids:

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)

        # Get output folder
        path_cwave_output = os.path.join(input_dir, probe_folder, 'kilosort2', 'cwaves')

        # Get sampling rate
        metafile_name = '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id)
        apbin_metafile_path = os.path.join(input_dir, probe_folder, metafile_name)
        ap_meta_dict = readSGLX.readMeta(pathlib.Path(apbin_metafile_path))
        imSampRate = float(ap_meta_dict['imSampRate'])  # probe-specific

        # Load mean waveform data from C_waves
        path_mean_waveforms = os.path.join(path_cwave_output, 'mean_waveforms.npy')
        if not os.path.isfile(path_mean_waveforms):
            print('Skipping probe. No mean waveforms at', path_mean_waveforms)
            continue
        mean_waveforms = np.load(path_mean_waveforms)

        # Get peak channels information
        #kilosort_folder = [f for f in os.listdir(os.path.join(input_dir, probe_folder)) if 'kilosort' in f][0]
        clus_info = pd.read_csv(os.path.join(input_dir, probe_folder, 'kilosort2', 'cluster_info.tsv'), sep='\\t')

        peak_channels = clus_info['ch'].values

        # Iterate over all clusters
        waveform_metrics_df = []
        for cluster_id in range(mean_waveforms.shape[0]):
            peak_chan = peak_channels[cluster_id]

            # Metrics for a single cluster waveform
            metrics_df, trough_idx, peak_idx = calculate_waveform_metrics_from_avg(
                avg_waveform=np.array(mean_waveforms[cluster_id, peak_chan, :]),
                cluster_id=cluster_id,
                peak_channel=peak_chan,
                sample_rate=imSampRate,
                upsampling_factor=250)
            metrics_df['peak_channel'] = peak_chan
            metrics_df['trough_idx'] = trough_idx #index of waveform trough
            metrics_df['peak_idx'] = peak_idx #index of waveform peak

            waveform_metrics_df.append(metrics_df)

        # Concatenate all cluster metrics into a single dataframe
        waveform_metrics_df = pd.concat(waveform_metrics_df)

        # Save dataframe
        waveform_metrics_df.to_csv(os.path.join(path_cwave_output, 'waveform_metrics.csv'), index=False)

    return






