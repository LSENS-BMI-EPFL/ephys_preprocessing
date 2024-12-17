#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: run_lfp_analysis.py
@time: 2/15/2024 11:13 PM
@description: Local field potential analysis for depth estimation and probe localization.
"""

# Imports
import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from loguru import logger
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from scipy.signal import welch
from scipy.ndimage.filters import gaussian_filter1d


# Import modules
from utils import readSGLX
from utils.sglx_meta_to_coords import MetaToCoords


def find_range(x, a, b, option='within'):
    """
    Find indices of data within or outside range [a,b].
    :param x:   numpy.ndarray  Data to search
    :param a:  float or int   Minimum value
    :param b:  float or int   Maximum value
    :param option:  String   'within' or 'outside'
    :return:    numpy.ndarray   Indices of x that fall within or outside specified range
    """
    if option == 'within':
        return np.where(np.logical_and(x >= a, x <= b))[0]
    elif option == 'outside':
        return np.where(np.logical_or(x < a, x > b))[0]
    else:
        raise ValueError('unrecognized option parameter: {}'.format(option))


def rms(data):
    """
    Computes root-mean-squared voltage of a signal.
    :param data: numpy.ndarray
    :return:  float  RMS value
    """
    return np.power(np.mean(np.power(data.astype('float32'), 2)), 0.5)


def printProgressBar(iteration, total, prefix='', suffix='', decimals=0, length=40, fill='*'):
    """
    Call in a loop to create terminal progress bar.
    Code from https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console.
    :param iteration:  int  Current iteration
    :param total:   int  Total iterations
    :param prefix:  str  (optional)  Prefix string
    :param suffix: str  (optional)  Suffix string
    :param decimals: int  (optional)  Positive number of decimals in percent complete
    :param length: int  (optional)  Character length of bar
    :param fill:  str  (optional)  Bar fill character
    :return:
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '░' * (length - filledLength)
    sys.stdout.write('\r%s %s %s%% %s' % (prefix, bar, percent, suffix))
    sys.stdout.flush()

    if iteration == total:
        print()


def find_surface_channel(lfp_data, ephys_params, params, xCoord, yCoord, shankInd):
    """
    Computes surface channel from LFP band data.
    Updated to use the site positions and estimate surface y (from tip) for shank 0.
    From: https://github.com/jenniferColonell/ecephys_spike_sorting/blob/master/ecephys_spike_sorting/modules/depth_estimation/depth_estimation.py#L79
    :param lfp_data: numpy.ndarray  (N samples x M channels)
    :param ephys_params:  dict  Ephys parameters
    :param params: dict  Processing parameters
    :param xCoord: numpy.ndarray  x-coordinates of probe sites
    :param yCoord: numpy.ndarray  y-coordinates of probe sites
    :param shankInd: numpy.ndarray  Shank index of probe sites
    :return: dict  Output dictionary with channel at brain surface and agar/air channels (approximate))
    """

    nchannels = ephys_params['num_channels']
    sample_frequency = ephys_params['lfp_sample_rate']

    # Get LFP data dimensions
    lfp_samples, lfp_channels = lfp_data.shape

    # Get parameters
    smoothing_amount = params['smoothing_amount']
    power_thresh = params['power_thresh']
    diff_thresh = params['diff_thresh']
    freq_range = params['freq_range_gamma'] # frequency range for gamma band
    freq_range_profile = params['freq_range_spiking'] # frequency range for LFP spiking profile
    saline_range = params['saline_range_um']
    nfft = params['nfft']
    n_passes = params['n_passes']

    # Initialize candidate lcoations
    candidates = np.zeros((n_passes,))

    # Determine number of samples per pass
    samples_per_pass = int(sample_frequency * (params['skip_s_per_pass'] + 1))
    max_passes = int(np.floor(lfp_samples / samples_per_pass))
    passes_used = min(n_passes, max_passes)

    # Use channels only on shank0, to yield a single estimate for the surface z
    channels = np.squeeze(np.asarray(np.where(shankInd == 0)))

    # Remove reference channels
    channels = np.delete(channels, ephys_params['reference_channels'])
    nchannels_used = channels.size

    # Get channels in saline
    chan_y = np.squeeze(yCoord[channels])
    in_saline_range = np.squeeze((chan_y > saline_range[0]) & (chan_y < saline_range[1]))
    saline_chan = np.where(in_saline_range)

    max_y = np.max(chan_y) # default init

    # Loop over passes
    for p in range(passes_used):

        startPt = int(sample_frequency * params['skip_s_per_pass'] * p)
        endPt = startPt + int(sample_frequency)

        chunk = np.copy(lfp_data[startPt:endPt, channels])

        # Subtract DC offset for all channels
        for ch in np.arange(nchannels_used):
            chunk[:, ch] = chunk[:, ch] - np.median(chunk[:, ch])

        # Reduce noise by correcting each timepoint with the signal in saline
        for ch in np.arange(nchannels_used):
            saline_chunk = np.squeeze(chunk[:, saline_chan])
            saline_median = np.median(saline_chunk, 1)
            chunk[:, ch] = chunk[:, ch] - saline_median

        # Init. power spectrum array
        power = np.zeros((int(nfft / 2 + 1), nchannels_used))

        # Compute power spectral density for each channel
        for ch in np.arange(nchannels_used):
            printProgressBar(p * nchannels_used + ch + 1, nchannels_used * n_passes)

            sample_frequencies, Pxx_den = welch(chunk[:, ch], fs=sample_frequency, nfft=nfft) #Welch method
            power[:, ch] = Pxx_den

        # Find indices of data within frequency range
        in_range = find_range(sample_frequencies, 0, params['max_freq'])
        in_range_gamma = find_range(sample_frequencies, freq_range[0], freq_range[1])
        in_range_spiking = find_range(sample_frequencies, freq_range_profile[0], freq_range_profile[1])

        # Compute mean of the log-power over input range
        values_gamma = np.log10(np.mean(power[in_range_gamma, :], 0))
        values_gamma = gaussian_filter1d(values_gamma, smoothing_amount)
        values_spiking = np.log10(np.mean(power[in_range_spiking, :], 0))
        values_spiking = gaussian_filter1d(values_spiking, smoothing_amount)

        # Find surface channels using power and derivative thresholds
        surface_channels = np.where((np.diff(values_gamma) < diff_thresh) * (values_gamma[:-1] < power_thresh))[0] # can return a list of channels
        surface_y = chan_y[surface_channels]

        # If no channels meet the criteria, use the maximum y position
        if len(surface_y > 0):
            candidates[p] = np.max(surface_y) # keep max y position if list given
        else:
            candidates[p] = max_y # else, take max y position of probe if no channels meet criteria

    # Determine brain surface channel taking median over iterations
    surface_y = np.median(candidates)

    # Determine air surface channel (separating agar/saline from air)
    air_y = np.min([surface_y + params['air_gap_um'], max_y]) # note: this is useless right now

    # Format and save results as json file
    output_dict = {
        'values_gamma': values_gamma.tolist(),
        'values_spiking':  values_spiking.tolist(),
        'surface_y': surface_y,
        'candidates': candidates.tolist(),
        'air_y': air_y
    }
    with open(os.path.join(params['figure_location'], 'depth.json'), 'w') as f:
        json.dump(output_dict, f)

    # Plot results
    if params['save_figure']:
        plot_results(chunk, power, in_range, values_gamma, nchannels_used, chan_y, surface_y, power_thresh, diff_thresh,
                     params['figure_location'])

    return output_dict

def plot_results(lfp_data_chunk,power,in_range,values,nchannels,chan_y,surface_y,power_thresh,diff_thresh,figure_location):
    """
    Plot results of LFP analysis.
    :param lfp_data_chunk: numpy.ndarray LFP data chunk
    :param power: numpy.ndarray Power spectrum
    :param in_range: numpy.ndarray Indices of data within frequency range for gamma
    :param values: numpy.ndarray Mean log-power over input range
    :param nchannels: int Number of channels
    :param chan_y: numpy.ndarray Y-coordinates of probe sites
    :param surface_y: float Y-coordinate of surface channel
    :param power_thresh: float Power threshold
    :param diff_thresh: float Power derivative threshold
    :param figure_location: str Path to figure location
    :return:
    """

    figsize = (5, 5)
    fig, axs = plt.subplots(2, 1, figsize=figsize, dpi=300)

    # Plot raw voltage data from last pass as image
    axs[0].imshow(np.flipud((lfp_data_chunk).T), aspect='auto')
    chunk_order = np.argsort(chan_y) # sort chunks by y position
    lfp_data_chunk[:, :] = lfp_data_chunk[:, chunk_order]
    # axs[0].imshow((chunk).T, aspect='auto',vmin=-1000,vmax=1000)
    axs[0].set_title('Raw voltage')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Channel number')
    axs[0].set_yticks(ticks=[384,0], labels=[1, 384])
    #axs[0].axhline(y=int((surface_y+175)/20), color='white', linestyle=':')

    # Plot the log(power) from lass pass as image
    axs[1].imshow(np.flipud(np.log10(power[in_range, :]).T), aspect='auto')
    power[:, :] = power[:, chunk_order] # sort chunks by y position
    # axs[1].imshow(np.log10(power[in_range,:]).T, aspect='auto')
    axs[1].set_title('Log10 power')
    axs[1].set_xlabel('Frequency (Hz)')
    axs[1].set_ylabel('Channel number')
    axs[1].set_yticks(ticks=[384,0], labels=[1, 384])
    #axs[1].axhline(y=int((surface_y+175)/20), color='white', linestyle=':')

    # Save figure
    fig.tight_layout()
    file_path = os.path.join(figure_location, 'probe_depth_matrix.png')
    fig.savefig(file_path, dpi='figure', bbox_inches='tight')

    figsize = (5, 5)
    fig, axs = plt.subplots(1, 2, figsize=figsize, dpi=300)

    # Plot power used in threshold calculation (=mean power over input range) vs. y-position
    y_sorted = chan_y[chunk_order]
    axs[0].plot(values[chunk_order], y_sorted)
    axs[0].plot([power_thresh, power_thresh], [chan_y[0], chan_y[nchannels - 1]], ls='--', c='k', label = 'power thresh.')
    #axs[0].axvline(x=surface_y, color='r', linestyle='--')
    axs[0].set_title('Mean power over input range')
    axs[0].set_ylabel(r'y-position from tip ($\mu$m)')
    axs[0].set_xlabel('Log10 power')

    try:
        surface_index = np.min(np.where(y_sorted > surface_y))
        axs[0].plot([-2, 2], [y_sorted[surface_index], y_sorted[surface_index]], ls='--', c='r', label = 'surface')
    except ValueError:
        pass

    axs[0].legend(frameon=False, loc='upper right')

    # Plot power derivative in threshold calculation vs. y-position

    #axs[1].plot(np.diff(values[chunk_order]), y_sorted[0:nchannels - 1])
    axs[1].plot(np.diff(values[chunk_order]), y_sorted[0:nchannels - 1]) # smoothed for visualization only
    axs[1].plot([diff_thresh, diff_thresh], [chan_y[0], chan_y[nchannels - 2]], ls='--', c='k', label = 'diff. thresh.')
    axs[1].plot([-0.2, diff_thresh], [surface_y, surface_y], ls='--', c='r', label = 'surface')
    axs[1].set_title('Surface channel: {}'.format(surface_y))
    axs[1].set_ylabel(r'y-position from tip ($\mu$m)')
    axs[1].set_xlabel('Log10 power derivative')
    axs[1].legend(frameon=False, loc='upper right')

    # Save figure
    fig.tight_layout()
    file_path = os.path.join(figure_location, 'probe_depth_profile.png')
    fig.savefig(file_path, dpi='figure', bbox_inches='tight')

    return

def bandpower(x, fs, fmin, fmax):
    f, Pxx = sp.signal.periodogram(x, fs=fs)
    ind_min = sp.argmax(f > fmin) - 1
    ind_max = sp.argmax(f > fmax) - 1
    return sp.trapz(Pxx[ind_min: ind_max], f[ind_min: ind_max])

def calculate_average_power(lfp_data, sampling_rate, freq_band):
    """
    Calculate average power in a specified frequency band.

    Parameters:
        lfp_data (ndarray): Array representing LFP data with dimensions (time, channels).
        sampling_rate (float): Sampling rate of the LFP data.
        freq_band (tuple): Frequency band of interest (start_freq, end_freq).

    Returns:
        float: Average power in the specified frequency band.
    """
    # Perform FFT along the time axis
    fft_result = np.fft.fft(lfp_data, axis=0)

    # Calculate corresponding frequencies
    freqs = np.fft.fftfreq(lfp_data.shape[0], 1 / sampling_rate)

    # Find indices corresponding to the specified frequency band
    start_index = np.argmax(freqs >= freq_band[0])
    end_index = np.argmax(freqs >= freq_band[1])

    # Calculate power spectrum
    power_spectrum = np.abs(fft_result) ** 2

    # Average power within the frequency band
    average_power = np.mean(np.sum(power_spectrum[start_index:end_index, :], axis=0))

    return average_power


def get_lfp_profile(lfp_data, ephys_params, params): #todo: not used, delete at some point?

    # Remove reference channels
    #lfp_data = np.delete(lfp_data, ephys_params['reference_channels'], axis=1)
    lfp_data = lfp_data[:, ::10]

    # Downsample
    lfp_data = lfp_data[::int(ephys_params['lfp_sample_rate']/1000), :]
    #lfp_data = lfp_data.mean(axis=1)

    freq_low = 500
    freq_high = 1250
    # Calculate average power in the frequency range of interest
    #lfp_band_power = bandpower(lfp_data, 1, freq_low, freq_high)
    fs = ephys_params['lfp_sample_rate']
    lfp_band_power = calculate_average_power(lfp_data=lfp_data,
                                             sampling_rate=ephys_params['lfp_sample_rate']/1000,
                                             freq_band=(freq_low, freq_high))

    # Smooth depth profile across number of channels
    n_chan_smooth = 10
    lfp_band_power_smooth = np.convolve(lfp_band_power, np.ones(n_chan_smooth) / n_chan_smooth, mode='same')

    # Plot to see
    #plt.figure()
    #plt.plot(lfp_band_power)
    #plt.plot(lfp_band_power_smooth)
    #plt.show()


    return


def main(input_dir):
    """
    Analyze local field potential data to estimate probe depth and localization.
    :param input_dir: path to CatGT preprocessed data
    :return:

    """

    # Processing parameters
    params = {
        "hi_noise_thresh": 50.0, #default is 50.0
        "lo_noise_thresh": 3.0, #default is 3.0
        "save_figure": True,
        "smoothing_amount": 5, #default is 5
        "power_thresh": 2.5, #default is 2.5
        "diff_thresh": -0.06, #default is -0.06
        "freq_range_gamma": [0, 10], #default is [0, 10] #Todo use real gamma?
        "freq_range_spiking": [500, 1250],
        "max_freq": 150, #default is 150
        "saline_range_um": [2000, 4000], # overwritten using insertion metadata
        "n_passes": 10, #default is 10
        "air_gap_um": 1000, #default is 1000
        "time_interval": 5, #default is 5
        "skip_s_per_pass": 10,  #default is 10
        "start_time": 10, #default is 5
        "nfft": 4096  # default but why
    }

    # Get epoch name and probe folders
    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Path to probe insertion info
    probe_info_path = r'M:\analysis\Axel_Bisi\mice_info\probe_insertion_info.xlsx'
    probe_info_df = pd.read_excel(probe_info_path)

    # Perform computations for each probe separately
    for probe_id, probe_folder in zip(probe_ids, probe_folders):

        # Get probe insertion depth for saline range
        mouse_name = epoch_name.split('_')[0]
        probe_row = probe_info_df[(probe_info_df['mouse_name'] == mouse_name)
                                  &
                                  (probe_info_df['probe_id'] == int(probe_id))]
        insertion_depth = probe_row['depth'].values[0]

        # Check that depth estimation using LFP is possible using electrodes in saline
        if insertion_depth > 4000:
            logger.warning('Probe {} insertion depth ({}) is too deep for LFP depth estimation. Skipped.'.format(probe_id, insertion_depth))
            continue

        params['saline_range_um'] = [insertion_depth, 4000]

        # Path to probe and output
        path_to_probe_folder = os.path.join(input_dir, probe_folder)
        path_to_output = os.path.join(path_to_probe_folder, 'depth')
        if not os.path.exists(path_to_output):
            os.makedirs(path_to_output)
        params['figure_location'] = os.path.join(path_to_output)

        # Path to metadata files
        ap_meta_file = os.path.join(input_dir, probe_folder,
                                    '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id))  # for channel coordinates
        lfp_meta_file = os.path.join(path_to_probe_folder,
                                     '{}_tcat.imec{}.lf.meta'.format(epoch_name, probe_id))  # for LFP info

        # Get LFP metadata
        lfp_meta = readSGLX.readMeta(Path(lfp_meta_file))
        num_channels = int(lfp_meta['nSavedChans'])
        lfp_sample_rate = float(lfp_meta['imSampRate'])
        reference_channels = [191]  # default for Neuropixels 1.0 probes
        ephys_params = {
            'num_channels': num_channels,
            'lfp_sample_rate': lfp_sample_rate,
            'reference_channels': reference_channels
        }

        # Load LFP band data and reshape
        lfp_data_file = os.path.join(path_to_probe_folder, '{}_tcat.imec{}.lf.bin'.format(epoch_name, probe_id))
        raw_data_lfp = np.memmap(lfp_data_file, dtype='int16', mode='r')
        lfp_data = np.reshape(raw_data_lfp, (int(raw_data_lfp.size / num_channels), num_channels))

        # Get probe channel coordinates
        coords = MetaToCoords(metaFullPath=Path(ap_meta_file),
                              outType=-1,
                              badChan=np.zeros((0), dtype='int'),
                              destFullPath='',
                              showPlot=False)
        xCoord = coords[0]
        yCoord = coords[1]
        shankInd = coords[2]

        # -----------------------
        # Compute surface channel
        # -----------------------

        logger.info('Computing surface channel for IMEC probe {}.'.format(probe_id))
        output_lfp = find_surface_channel(lfp_data=lfp_data, ephys_params=ephys_params, params=params, xCoord=xCoord, yCoord=yCoord, shankInd=shankInd)

        # Plot  values in spiking range profile, starting from the surface channel
        #print('Plotting LFP profile...')
        lfp_values = np.array(output_lfp['values_spiking'])
        #lfp_values_to_plot = lfp_values[int(output_lfp['surface_y']/20):-1]
        surface_y = output_lfp['surface_y']
        #print(lfp_values.shape, surface_y)
        fig, ax = plt.subplots(1,1, figsize=(5,10), dpi=300)
        x_range = np.arange(min(lfp_values), max(lfp_values), lfp_values.shape)
        x_range = np.linspace(min(lfp_values), max(lfp_values), lfp_values.shape[0])
        #print(x_range.shape, lfp_values[int(surface_y/20):-1].shape)
        #ax.plot(lfp_values[int(surface_y/20):-1], x_range, lw=2, c='k')
        ax.plot(lfp_values, x_range, lw=2, c='k')
        ax.axhline(y=surface_y/20, color='r', linestyle='--')
        ax.set_title('imec{}'.format(probe_id))
        ax.set_xlabel('Mean log10(power)')
        ax.set_ylabel('Channel number')
        #plt.show()

        # -------------------
        # Compute LFP profile
        # -------------------
        #print('Computing LFP profile...')
        #lfp_profile = get_lfp_profile(lfp_data=lfp_data, ephys_params=ephys_params, params=params)






    return
