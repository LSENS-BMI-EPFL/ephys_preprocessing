#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: ephys_datajoint
@file: waveform_metrics_utils.py
@time: 28.04.2023 09:10
@description: Waveform metrics calculation.
"""

# Imports
import numpy as np
import pandas as pd

from scipy.stats import linregress
from scipy.signal import resample


def calculate_waveform_metrics_from_avg(avg_waveform,
                                        cluster_id,
                                        peak_channel,
                                        sample_rate,
                                        upsampling_factor=200):
    """
    Calculate metrics for an array of waveforms for a single cluster.

    Metrics come from Jia et al. (2019) High-density extracellular probes reveal
    dendritic backpropagation and facilitate neuron classification. J Neurophys
    https://doi.org/10.1152/jn.00680.2018


    Inputs:
    -------
    avg_waveform : numpy.ndarray (num_channels x num_samples)
        from C_waves output
    cluster_id : int
        ID for cluster
    peak_channel : int
        Location of waveform peak, read from C_waves input clus_Table
    channel_map : numpy.ndarray
        Channels used for spike sorting
    sample_rate : float
        Sample rate in Hz
    upsampling_factor : float
        Relative rate at which to upsample the spike waveform
    spread_threshold : float
        Threshold for computing spread of 2D waveform
    site_range : float
        Number of sites to use for 2D waveform metrics
    site_x, site_y : channel positions in um

    Outputs:
    -------
    metrics : pandas.DataFrame
        Single-row table containing all metrics

    """

    # Upsample mean waveform
    num_samples = avg_waveform.shape[0]
    new_sample_count = int(num_samples * upsampling_factor)
    avg_waveform = resample(avg_waveform, new_sample_count)

    timestamps = np.linspace(0, num_samples / sample_rate, new_sample_count)

    # Calculate metrics
    duration, trough_idx, peak_idx = calculate_waveform_duration(avg_waveform, timestamps) #trough-to-peak
    halfwidth = calculate_waveform_halfwidth(avg_waveform, timestamps)
    pt_ratio = calculate_waveform_PT_ratio(avg_waveform)
    repolarization_slope = calculate_waveform_repolarization_slope(
        avg_waveform, timestamps)
    recovery_slope = calculate_waveform_recovery_slope(
        avg_waveform, timestamps)

    # Format as one-row dataframe
    data = [[cluster_id, peak_channel, duration, halfwidth, pt_ratio, repolarization_slope,
             recovery_slope]]

    metrics = pd.DataFrame(data, columns=['cluster_id', 'peak_channel', 'duration', 'halfwidth',
                                    'pt_ratio', 'repolarization_slope', 'recovery_slope'])

    return metrics, trough_idx, peak_idx

# ==========================================================

# MEAN WAVEFORM HELPERS

# ==========================================================


def calculate_waveform_duration(waveform, timestamps):
    """
    Duration (in seconds) between peak and trough

    Inputs:
    ------
    waveform : numpy.ndarray (N samples)
    timestamps : numpy.ndarray (N samples)

    Outputs:
    --------
    duration : waveform duration in milliseconds

    """

    trough_idx = np.argmin(waveform)
    peak_idx = np.argmax(waveform)
    if peak_idx < trough_idx:
        peak_idx = trough_idx + np.argmax(waveform[trough_idx:-1]) # get peak index, starting from trough


    ## to avoid detecting peak before trough
    #if waveform[peak_idx] > np.abs(waveform[trough_idx]):
    #    duration = timestamps[peak_idx:][np.where(waveform[peak_idx:] == np.min(waveform[peak_idx:]))[0][0]] - \
    #               timestamps[peak_idx]
    #    #new_peak_idx = np.where(waveform[peak_idx:] == np.min(waveform[peak_idx:]))[0][0]
    #else:
    #    duration = timestamps[trough_idx:][np.where(waveform[trough_idx:] == np.max(waveform[trough_idx:]))[0][0]] - \
    #               timestamps[trough_idx]
    #
    try:
        duration = timestamps[peak_idx] - timestamps[trough_idx]
    except IndexError:
        duration = np.nan

    if duration > 1.2:
        duration = np.nan

    return duration * 1e3, trough_idx, peak_idx


def calculate_waveform_halfwidth(waveform, timestamps):
    """
    Spike width (in seconds) at half max amplitude

    Inputs:
    ------
    waveform : numpy.ndarray (N samples)
    timestamps : numpy.ndarray (N samples)

    Outputs:
    --------
    halfwidth : waveform halfwidth in milliseconds

    """

    trough_idx = np.argmin(waveform)
    peak_idx = np.argmax(waveform)

    try:
        if waveform[peak_idx] > np.abs(waveform[trough_idx]):
            threshold = waveform[peak_idx] * 0.5
            thresh_crossing_1 = np.min(
                np.where(waveform[:peak_idx] > threshold)[0])
            thresh_crossing_2 = np.min(
                np.where(waveform[peak_idx:] < threshold)[0]) + peak_idx
        else:
            threshold = waveform[trough_idx] * 0.5
            thresh_crossing_1 = np.min(
                np.where(waveform[:trough_idx] < threshold)[0])
            thresh_crossing_2 = np.min(
                np.where(waveform[trough_idx:] > threshold)[0]) + trough_idx

        halfwidth = (timestamps[thresh_crossing_2] - timestamps[thresh_crossing_1])

    except ValueError as err:
        print(err)
        halfwidth = np.nan

    return halfwidth * 1e3


def calculate_waveform_PT_ratio(waveform):
    """
    Peak-to-trough ratio of 1D waveform

    Inputs:
    ------
    waveform : numpy.ndarray (N samples)

    Outputs:
    --------
    PT_ratio : waveform peak-to-trough ratio

    """

    trough_idx = np.argmin(waveform)

    peak_idx = np.argmax(waveform)

    PT_ratio = np.abs(waveform[peak_idx] / waveform[trough_idx])

    return PT_ratio


def calculate_waveform_repolarization_slope(waveform, timestamps, window=20):
    """
    Spike repolarization slope (after maximum deflection point)

    Inputs:
    ------
    waveform : numpy.ndarray (N samples)
    timestamps : numpy.ndarray (N samples)
    window : int
        Window (in samples) for linear regression

    Outputs:
    --------
    repolarization_slope : slope of return to baseline (V / s)

    """

    max_point = np.argmax(np.abs(waveform))

    waveform = - waveform * (np.sign(waveform[max_point]))  # invert if we're using the peak

    repolarization_slope = linregress(timestamps[max_point:max_point + window], waveform[max_point:max_point + window])[
        0]

    return repolarization_slope * 1e-6


def calculate_waveform_recovery_slope(waveform, timestamps, window=20):
    """
    Spike recovery slope (after repolarization)

    Inputs:
    ------
    waveform : numpy.ndarray (N samples)
    timestamps : numpy.ndarray (N samples)
    window : int
        Window (in samples) for linear regression

    Outputs:
    --------
    recovery_slope : slope of recovery period (V / s)

    """

    max_point = np.argmax(np.abs(waveform))

    waveform = - waveform * (np.sign(waveform[max_point]))  # invert if we're using the peak

    peak_idx = np.argmax(waveform[max_point:]) + max_point

    recovery_slope = linregress(timestamps[peak_idx:peak_idx + window], waveform[peak_idx:peak_idx + window])[0]

    return recovery_slope * 1e-6
