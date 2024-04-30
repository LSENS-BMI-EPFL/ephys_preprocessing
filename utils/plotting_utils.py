#! /usr/bin/env/python3
"""
@author: Axel Bisi
@project: EphysUtils
@file: plotting_utils.py
@time: 10.05.2023 10:47
"""

# Imports
import numpy as np
from matplotlib import colors
import matplotlib.colors as mc
import colorsys
import scipy.ndimage


def remove_top_right_frame(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return

def color_to_rgba(color_name):
    """
    Converts color name to RGB.
    :param color_name:
    :return:
    """

    return colors.to_rgba(color_name)

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    From: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    :param color: Matplotlib color string.
    :param amount: Number between 0 and 1.
    :return:
    """

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def adjust_lightness(color, amount=0.5):
    """
    Same as lighten_color but adjusts brightness to lighter color if amount>1 or darker if amount<1.
    Input can be matplotlib color string, hex string, or RGB tuple.
    From: https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    :param color: Matplotlib color string.
    :param amount: Number between 0 and 1.
    :return:
    """

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def make_cmap_n_from_color_lite2dark(color, N):
    """
    Make ListedColormap from matplotlib color of size N using the lighten_color function.
    :param color: Matplotlib color string.
    :param N: Number of colors to have in cmap.
    :return:
    """
    light_factors = np.linspace(0.2, 1, N)
    cmap = colors.ListedColormap(colors=[lighten_color(color, amount=i) for i in light_factors])
    return cmap

def halfgaussian_kernel1d(sigma, radius):
    """
    Computes a 1-D Half-Gaussian convolution kernel.
    From: https://stackoverflow.com/questions/71003634/applying-a-half-gaussian-filter-to-binned-time-series-data-in-python
    """
    sigma2 = sigma * sigma
    x = np.arange(0, radius+1)
    phi_x = np.exp(-0.5 / sigma2 * x ** 2)
    phi_x = phi_x / phi_x.sum()

    return phi_x

def halfgaussian_filter1d(input, sigma, axis=-1, output=None,
                      mode="constant", cval=0.0, truncate=4.0):
    """
    Convolves a 1-D Half-Gaussian convolution kernel.
    From: https://stackoverflow.com/questions/71003634/applying-a-half-gaussian-filter-to-binned-time-series-data-in-python
    """
    sd = float(sigma)
    # make the radius of the filter equal to truncate standard deviations
    lw = int(truncate * sd + 0.5)
    weights = halfgaussian_kernel1d(sigma, lw)
    origin = -lw // 2
    return scipy.ndimage.convolve1d(input, weights, axis, output, mode, cval, origin)

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