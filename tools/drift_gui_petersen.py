from qtpy import QtGui
from qtpy.QtWidgets import QWidget, QApplication, QSizePolicy, QMainWindow, QGridLayout
from qtpy.QtCore import Qt
from qtpy.QtGui import QTransform

import matplotlib

import pyqtgraph as pg
import numpy as np

from iblutil.numerical import bincount2D
import one.alf.io as alfio

import matplotlib as mpl
import pandas as pd
from pathlib import Path

import spikeinterface.full as si
import spikeinterface.preprocessing as sp
from spikeinterface.postprocessing.unit_locations import dtype_localize_by_method
from spikeinterface.extractors.alfsortingextractor import ALFSortingExtractor
from spikeinterface.metrics.quality import ComputeQualityMetrics
from spikeinterface.sortingcomponents.tools import (
    create_sorting_analyzer_with_existing_templates,
)

#%%
# CHANGE PATHS AS NECESSARY
RAW_PATH = Path(r"M:\analysis\Axel_Bisi\data\AB163\AB163_20250419_155630\Ephys\catgt_AB163_g0\AB163_g0_imec0")
STREAM_NAME = 'imec0.ap'
SORTING_PATH = Path(r"M:\analysis\Axel_Bisi\data\AB163\AB163_20250419_155630\Ephys\catgt_AB163_g0\AB163_g0_imec0\ibl_format")
METRICS_FILE = Path(r"M:\analysis\Axel_Bisi\data\AB163\AB163_20250419_155630\Ephys\catgt_AB163_g0\AB163_g0_imec0\kilosort2\cluster_bc_unitType.tsv")
MOTION_PATH = Path(r'M:\analysis\Axel_Bisi\data\AB163\AB163_20250419_155630\Ephys\catgt_AB163_g0\AB163_g0_imec0\dredge_fast\motion')

# DRIFT INTERVALS DURING RECORDING TO COMPARE
INTERVALS = [
    [0, 2000],
    [2000, 4000],
]

#%%
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


GOOD_RGB = np.array([40, 170, 70, 220], dtype=np.uint8)
MUA_RGB = np.array([0, 0, 255, 220], dtype=np.uint8)
BAD_RGB = np.array([210, 50, 50, 220], dtype=np.uint8)
NOLABEL_RGB = np.array([70, 70, 70, 220], dtype=np.uint8)

FONT = QtGui.QFont()
FONT.setPointSize(8)

DEFAULT_PEN = pg.mkPen(None)
SELECTED_PEN = pg.mkPen(255, 220, 0, width=6)
DEFAULT_SIZE = 10
SELECTED_SIZE = 15

INTERVAL_COLORS = [
    (230, 159, 0, 70),
    (86, 180, 233, 70),
    (0, 158, 115, 70),
]


#%%
def load_data_with_spike_interface():

    # Load in the raw data
    recording = si.read_spikeglx(folder_path=RAW_PATH, stream_name=STREAM_NAME)
    recording.reset_times()
    fs = recording.sampling_frequency
    sorting = ALFSortingExtractor(SORTING_PATH, fs)

    # Load in the spike sorting data
    spikes = alfio.load_object(SORTING_PATH, "spikes")
    templates = alfio.load_object(SORTING_PATH, "templates")
    templates_array = np.nan_to_num(templates["waveforms"])
    winds = templates["waveformsChannels"]

    sparsity = np.zeros((len(templates_array), 384), dtype=bool)
    for i in range(len(templates_array)):
        mask = winds[i] < 384
        sparsity[i, winds[i, mask]] = True

    templates = si.Templates(
        templates_array * 1e6,
        sampling_frequency=fs,
        nbefore=templates_array.shape[2] // 2,
        is_in_uV=True,
        sparsity_mask=sparsity,
        probe=recording.get_probe(),
        channel_ids=recording.channel_ids,
    )

    positions = np.zeros(
        len(spikes["depths"]), dtype=dtype_localize_by_method["center_of_mass"]
    )
    positions["y"] = spikes["depths"]


    analyzer = create_sorting_analyzer_with_existing_templates(
        sorting,
        recording,
        spike_amplitudes=1e6 * spikes["amps"],
        spike_locations=positions,
        noise_levels=np.ones(recording.get_num_channels()),
        templates=templates,
    )

    # Load in the metrics
    metrics = pd.read_csv(METRICS_FILE, sep='\t')
    metrics['gui_label'] = metrics['bc_unitType']
    analyzer.extensions["quality_metrics"] = ComputeQualityMetrics(analyzer)
    analyzer.extensions["quality_metrics"].set_data('metrics', metrics)
    analyzer.extensions["quality_metrics"].run_info["run_completed"] = True
    analyzer.extensions["quality_metrics"].run_info["runtime_s"] = 0

    # Load in the motion data
    motion = sp.load_motion_info(MOTION_PATH)['motion']

    return analyzer, motion



class MainGUI(QMainWindow):

    # --------------------------------------------------------------------------------------------
    # Main
    # --------------------------------------------------------------------------------------------

    def __init__(self):
        super().__init__()

        self.analyzer, self.motion = load_data_with_spike_interface()

        if INTERVALS is None:
            self.intervals = [
                [0, self.analyzer.get_total_duration() //2],
                [self.analyzer.get_total_duration()//2, self.analyzer.get_total_duration()],
            ]
        else:
            self.intervals = INTERVALS

        # Initialize the setu[
        self.init_gui()

        # Load the data
        self.compute_average_metrics_by_cluster()
        has_label = 'gui_label' in self.analyzer.get_extension('quality_metrics').get_data()
        if not has_label:
            self.cluster_rgba = np.tile(NOLABEL_RGB, (self.clust_idx.size, 1))
        else:
            self.cluster_rgba = np.tile(BAD_RGB, (self.clust_idx.size, 1))
            mua_idx = self.analyzer.get_extension('quality_metrics').get_data()['gui_label'] == 'MUA'
            self.cluster_rgba[mua_idx] = MUA_RGB
            good_idx = self.analyzer.get_extension('quality_metrics').get_data()['gui_label'] == 'GOOD'
            self.cluster_rgba[good_idx] = GOOD_RGB

        # Initialise plots
        self.plot_depth_amp_scatter()
        self.plot_drift_overview()
        self.plot_spike_raster()
        self.plot_drift_lines(self.fig2, self.motion.spatial_bins_um)
        self.plot_drift_lines(self.fig3, np.arange(self.motion.spatial_bins_um.size) * 20)
        self.plot_drift_intervals()

        # Initialise cluster cycle mode
        self.spot_by_cluster_idx = {int(spot.data()): spot for spot in self.scatter.points()}
        self.cycle_mode = 'all'
        self.cycle_cluster_ids = self.clust_idx

        # Initialise the first selected cluster
        first_spot = self.scatter.points()[0]
        self.selected_clust_idx = int(first_spot.data())
        self.selected_spot = None

        self.set_selected_cluster(self.selected_clust_idx)

    # --------------------------------------------------------------------------------------------
    # Initialization and setup
    # --------------------------------------------------------------------------------------------
    def init_gui(self):
        self.resize(1800, 1000)
        self.setWindowTitle('Drift QC GUI')
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        central_widget = QWidget()
        central_layout = QGridLayout(central_widget)
        self.setCentralWidget(central_widget)
        self.central_widget = central_widget

        # Figure 1 scatter depth vs amp
        self.fig1 = pg.PlotWidget()
        self.fig1.setLabel('left', 'Depth along probe (um)')
        self.fig1.setLabel('bottom', 'Amplitude (uV)')
        self.fig1.getAxis('left').label.setFont(FONT)
        self.fig1.getAxis('bottom').label.setFont(FONT)

        # Figure 2 image of depth and firing rate
        self.fig2 = pg.PlotWidget()
        self.fig2.setLabel('left', 'Depth along probe (um)')
        self.fig2.setLabel('bottom', 'Time (s)')
        self.fig2.getAxis('left').label.setFont(FONT)
        self.fig2.getAxis('bottom').label.setFont(FONT)

        max_width = 80
        self.fig3 = pg.PlotItem()
        self.fig3.setTitle("<span style='font-size:8pt'>Drift overview</span>")
        self.fig3.getAxis('left').setWidth(max_width)
        self.fig3.setLabel('left', 'Drift (um)')
        self.fig3.getAxis('left').label.setFont(FONT)

        self.fig4 = pg.PlotItem()
        self.amp_img = pg.ImageItem()
        self.fig4.addItem(self.amp_img)
        self.fig4.setTitle("<span style='font-size:8pt'>Amplitude distribution across session</span>")
        self.fig4.getAxis('left').setWidth(max_width)
        self.fig4.setLabel('left', 'Amplitude (uV)', size=2)
        self.fig4.getAxis('left').label.setFont(FONT)

        self.fig5 = pg.PlotItem()
        self.fr_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2), name='fr')
        self.fig5.addItem(self.fr_curve)
        self.fig5.getAxis('left').setWidth(max_width)
        self.fig5.setTitle("<span style='font-size:8pt'>Mean cluster firing rate across session</span>")
        self.fig5.setLabel('left', 'Mean firing rate (spikes/s)')
        self.fig5.getAxis('left').label.setFont(FONT)

        self.fig6 = pg.PlotItem()
        self.amp_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2), name='amp')
        self.fig6.addItem(self.amp_curve)
        self.fig6.getAxis('left').setWidth(max_width)
        self.fig6.setTitle("<span style='font-size:8pt'>Mean amplitude across session</span>")
        self.fig6.setXRange(0, self.motion.temporal_bins_s[0].max() + 100)
        self.fig6.setLabel('left', 'Mean amplitude (uV)')
        self.fig6.getAxis('left').label.setFont(FONT)

        self.fig7 = pg.PlotItem()
        self.drift_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2), name='drift')
        self.fig7.addItem(self.drift_curve)
        self.fig7.getAxis('left').setWidth(max_width)
        self.fig7.setTitle("<span style='font-size:8pt'>Estimated drift at cluster depth across session</span>")
        self.fig7.setLabel('left', 'Estimated drift (um)')
        self.fig7.setLabel('bottom', 'Time (s)')
        self.fig7.getAxis('left').label.setFont(FONT)
        self.fig7.getAxis('bottom').label.setFont(FONT)

        self.fig8 = pg.PlotItem()
        self.fig8_hist_items = []
        self.fig8.getAxis('left').setWidth(max_width)
        self.fig8.setTitle("<span style='font-size:8pt'>Histogram of mean firing-rate per interval</span>")
        self.fig8.setLabel('left', 'Probability density')
        self.fig8.setLabel('bottom', 'Mean firing rate (spikes/bin)')
        self.fig8.getAxis('left').label.setFont(FONT)
        self.fig8.getAxis('bottom').label.setFont(FONT)

        self.fig9 = pg.PlotItem()
        self.fig9_hist_items = []
        self.fig9.getAxis('left').setWidth(max_width)
        self.fig9.setTitle("<span style='font-size:8pt'>Histogram of mean amplitude per interval</span>")
        self.fig9.setLabel('left', 'Probability density')
        self.fig9.setLabel('bottom', 'Mean amplitude (uV)')
        self.fig9.getAxis('left').label.setFont(FONT)
        self.fig9.getAxis('bottom').label.setFont(FONT)

        layout = pg.GraphicsLayout()
        widget = pg.GraphicsLayoutWidget()
        hist_layout = pg.GraphicsLayout()
        hist_widget = pg.GraphicsLayoutWidget()

        layout.addItem(self.fig3, row=0, col=0, rowspan=3)
        layout.addItem(self.fig4, row=4, col=0, rowspan=3)
        layout.addItem(self.fig5, row=7, col=0)
        layout.addItem(self.fig6, row=8, col=0)
        layout.addItem(self.fig7, row=9, col=0)

        hist_layout.addItem(self.fig8, row=0, col=0)
        hist_layout.addItem(self.fig9, row=1, col=0)

        widget.addItem(layout)
        hist_widget.addItem(hist_layout)

        central_layout.addWidget(self.fig1, 0, 0, 5, 1)
        central_layout.addWidget(self.fig2, 0, 1, 5, 1)
        central_layout.addWidget(widget, 0, 2, 5, 1)
        central_layout.addWidget(hist_widget, 1, 3, 3, 1)

    # --------------------------------------------------------------------------------------------
    # Load data
    # --------------------------------------------------------------------------------------------

    def compute_average_metrics_by_cluster(self):
        """Compute average amplitude, firing rate, and depth for each cluster."""
        df = pd.DataFrame()
        df['clusters'] = self.analyzer.sorting.to_spike_vector()["unit_index"]
        df['amps'] = self.analyzer.get_extension("spike_amplitudes").data['amplitudes']
        df['depths'] = self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"]

        avgs = df.groupby('clusters').agg(['mean', 'count'])
        self.clust_idx = avgs.index.values
        self.avg_amp = avgs['amps']['mean'].values
        self.avg_fr = avgs['depths']['count'].values / self.analyzer.get_total_duration()
        self.avg_depth = avgs['depths']['mean'].values

        del df, avgs

    def compute_data_for_cluster(self, clust_idx: int):

        idx = self.analyzer.sorting.to_spike_vector()["unit_index"] == clust_idx
        self.clust_times = self.analyzer.recording.get_times()[self.analyzer.sorting.to_spike_vector()["sample_index"][idx]]
        self.clust_amps = self.analyzer.get_extension("spike_amplitudes").data['amplitudes'][idx]


        # Get the mean amplitude and firing rate
        n_bins = self.motion.temporal_bins_s[0].size
        bin_indices = np.searchsorted(self.motion.temporal_bin_edges_s[0], self.clust_times, side="right") - 1
        bin_indices = np.clip(bin_indices, 0, n_bins - 1)

        amp_sum = np.bincount(bin_indices, weights=self.clust_amps, minlength=n_bins)
        amp_count = np.bincount(bin_indices, minlength=n_bins)

        self.mean_amp = np.full(n_bins, np.nan, dtype=float)
        valid = amp_count > 0
        self.mean_amp[valid] = amp_sum[valid] / amp_count[valid]
        self.mean_fr = np.full(n_bins, np.nan, dtype=float)
        self.mean_fr[valid] = amp_count[valid]

    # --------------------------------------------------------------------------------------------
    # Plots
    # --------------------------------------------------------------------------------------------
    def plot_depth_amp_scatter(self):

        self.scatter = pg.ScatterPlotItem()
        self.scatter.sigClicked.connect(self._on_plot_click)
        self.fig1.addItem(self.scatter)

        self.default_pen = pg.mkPen(None)
        self.selected_pen = pg.mkPen(255, 220, 0, width=6)
        self.default_size = 10
        self.selected_size = 15

        self.scatter.setData(
            x=self.avg_amp,
            y=self.avg_depth,
            data=self.clust_idx,
            brush=[pg.mkBrush(*rgba) for rgba in self.cluster_rgba],
            pen=DEFAULT_PEN,
            size=DEFAULT_SIZE,
            pxMode=True,
        )
        self.fig1.setXRange(-50, np.nanmin([np.nanmax(self.avg_amp) + 200, 1200]))
        self.fig1.setYRange(0, 4000)

    def plot_spike_raster(self):

        isnan = ~np.isnan(self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"])
        depths = self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"][isnan]
        times = self.analyzer.recording.get_times()[self.analyzer.sorting.to_spike_vector()["sample_index"][isnan]]

        fr, times, depths = bincount2D(times, depths, xbin=0.05, ybin=10, ylim=[np.min([0, np.min(depths)]), np.max([3840, np.max(depths)])])
        fr = fr.T

        xscale = (times[-1] - times[0]) / fr.shape[0]
        yscale = (depths[-1] - depths[0]) / fr.shape[1]
        yoffset = depths[0]
        xoffset = times[0]

        levels = np.quantile(np.mean(fr, axis=0), [0, 1])
        fr_img = pg.ImageItem()
        fr_img.setImage(fr)
        fr_img.setTransform(QTransform(xscale, 0.0, 0.0, 0.0, yscale, 0.0, xoffset, yoffset, 1.0))
        fr_img.setOpacity(0.8)
        cmap, lut, grad = get_color('binary')
        fr_img.setLookupTable(lut)
        fr_img.setLevels(levels)
        self.fig2.addItem(fr_img)
        self.fig2.setXRange(0, self.analyzer.get_total_duration() + 50)
        self.fig2.setYRange(0, 4000)

    def plot_drift_overview(self):
        # Compute drift interpolated at different depths
        drift_interp = self.motion.get_displacement_at_time_and_depth(
            times_s=self.motion.temporal_bins_s[0],
            locations_um=np.arange(0, 3840, 40),
            segment_index=0,
            grid=True,
        ).T

        xscale = (self.motion.temporal_bin_edges_s[0].max() - self.motion.temporal_bin_edges_s[0].min()) / self.motion.temporal_bins_s[0].size
        yscale = 3840 / drift_interp.shape[1]
        xoffset = self.motion.temporal_bin_edges_s[0].min()
        yoffset = 0

        drift_img = pg.ImageItem()
        drift_img.setImage(drift_interp)
        drift_img.setTransform(QTransform(xscale, 0.0, 0.0, 0.0, yscale, 0.0, xoffset, yoffset, 1.0))
        cmap, lut, grad = get_color('seismic')
        drift_img.setLookupTable(lut)
        drift_img.setLevels((-30, 30))

        self.fig2.addItem(drift_img)
        self.fig2.setXRange(0, self.analyzer.get_total_duration() + 50)
        self.fig2.setYRange(0, 4000)

    def plot_drift_lines(self, fig, offsets):

        for i, [disp, offset] in enumerate(zip(self.motion.displacement[0].T, offsets)):

            drift_plot = pg.PlotCurveItem()
            drift_plot.setData(x=self.motion.temporal_bins_s[0], y=disp + offset, pen='k', linewidth=4)
            fig.addItem(drift_plot)

        fig.setXRange(0, self.analyzer.get_total_duration() + 50)


    def plot_drift_intervals(self):

        t_min = float(np.nanmin(self.motion.temporal_bins_s[0]))
        t_max = float(np.nanmax(self.motion.temporal_bins_s[0]))
        for (start, end), color in zip(self.intervals, INTERVAL_COLORS):
            lo = max(float(start), t_min)
            hi = min(float(end), t_max)
            if hi <= lo:
                continue
            region = pg.LinearRegionItem(
                values=(lo, hi),
                orientation='vertical',
                movable=False,
                brush=pg.mkBrush(*color),
                pen=pg.mkPen(color[0], color[1], color[2], 130),
                swapMode='sort',
            )
            region.setZValue(-20)
            self.fig7.addItem(region)

    def plot_amplitude_distribution(self):

        self.amp_img.clear()
        bins = [400, 200]
        H, xedges, yedges = np.histogram2d(self.clust_times, self.clust_amps * 1e6, bins=bins)
        H[H == 0] = np.nan
        self.amp_img.setImage(H)
        xscale = (xedges[-1] - xedges[0]) / bins[0]
        yscale = (yedges[-1] - yedges[0]) / bins[1]
        xoffset = xedges[0]
        yoffset = yedges[0]
        self.amp_img.setTransform(QTransform(xscale, 0.0, 0.0, 0.0, yscale, 0.0, xoffset, yoffset, 1.0))
        cmap = pg.colormap.get("viridis")
        self.amp_img.setLookupTable(cmap.getLookupTable())
        self.fig4.setYRange(0, yedges[-1] + 20)
        self.fig4.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_mean_amplitude(self):

        self.amp_curve.clear()
        self.amp_curve.setData(self.motion.temporal_bins_s[0], self.mean_amp * 1e6)
        self.fig6.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_mean_firing_rate(self):
        self.fr_curve.clear()
        self.fr_curve.setData(self.motion.temporal_bins_s[0], self.mean_fr)
        self.fig5.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_estimated_depth(self):

        drift_at_depth = self.motion.get_displacement_at_time_and_depth(
            times_s=self.motion.temporal_bins_s[0],
            locations_um=[self.avg_depth[self.selected_clust_idx]],
            grid=True,
        )[0]

        self.drift_curve.clear()
        self.drift_curve.setData(self.motion.temporal_bins_s[0], drift_at_depth)
        self.fig7.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_histogram(self, data, fig, hist_items):

        for item in hist_items:
            fig.removeItem(item)

        hist_items = []

        data_max = float(np.nanmax(data))
        if data_max <= 0:
            data_max = 1.0
        bin_edges = np.linspace(0.0, data_max * 1.05, 24)
        bin_width = np.diff(bin_edges)
        centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        times = self.motion.temporal_bins_s[0]

        for (start, end), color in zip(self.intervals, INTERVAL_COLORS):
            interval_mask = (times >= start) & (times < end)
            vals = data[interval_mask]
            if vals.size == 0:
                continue

            counts, _ = np.histogram(vals, bins=bin_edges)
            total = counts.sum()
            if total == 0:
                continue
            counts = counts.astype(float) / (total * bin_width)
            bar = pg.BarGraphItem(
                x=centers,
                height=counts,
                width=bin_width * 0.92,
                brush=pg.mkBrush(color[0], color[1], color[2], 170),
                pen=pg.mkPen(color[0], color[1], color[2], 220),
            )
            fig.addItem(bar)
            hist_items.append(bar)

        return hist_items

    def plot_amplitude_histogram(self):

        self.fig9_hist_items = self.plot_histogram(self.mean_amp, self.fig9, self.fig9_hist_items)

    def plot_firing_rate_histogram(self):

        self.fig8_hist_items = self.plot_histogram(self.mean_fr, self.fig8, self.fig8_hist_items)


    # --------------------------------------------------------------------------------------------
    # Interactions
    # --------------------------------------------------------------------------------------------

    def _on_plot_click(self, _, point):
        clicked_pos = point[0].pos()
        clicked_spots = self.scatter.pointsAt(clicked_pos)
        if len(clicked_spots) == 0:
            return

        clicked_spot = clicked_spots[0]
        self.set_selected_cluster(int(clicked_spot.data()))

    def keyPressEvent(self, event):
        if event.modifiers() & Qt.ShiftModifier:
            mode_key_map = {
                Qt.Key_G: 'good',
                Qt.Key_M: 'mua',
                Qt.Key_B: 'bad',
                Qt.Key_A: 'all',
            }
            mode = mode_key_map.get(event.key())
            if mode is not None:
                self._set_cycle_mode(mode)
                event.accept()
                return

        if event.key() in (Qt.Key_Left, Qt.Key_Right):
            self._sync_cycle_ids_with_selected()
            if self.cycle_cluster_ids.size == 0:
                event.accept()
                return
            current_pos = int(np.where(self.cycle_cluster_ids == self.selected_clust_idx)[0][0])
            step = -1 if event.key() == Qt.Key_Left else 1
            next_pos = (current_pos + step) % self.cycle_cluster_ids.size
            self.set_selected_cluster(int(self.cycle_cluster_ids[next_pos]))
            event.accept()
            return

        super().keyPressEvent(event)

    # --------------------------------------------------------------------------------------------
    # Keyboard interaction for cycling through clusters based on quality labels
    # --------------------------------------------------------------------------------------------

    def _clusters_for_cycle_mode(self, mode: str) -> np.ndarray:
        """Return cluster IDs eligible for keyboard cycling under a mode."""
        if mode == 'good':
            return self.clust_idx[np.all(self.cluster_rgba == GOOD_RGB, axis=1)]
        if mode == 'mua':
            return self.clust_idx[np.all(self.cluster_rgba == MUA_RGB, axis=1)]
        if mode == 'bad':
            return self.clust_idx[np.all(self.cluster_rgba == BAD_RGB, axis=1)]

        return self.clust_idx

    def _set_cycle_mode(self, mode: str):
        """Set active keyboard cycle mode and keep selection valid for that mode."""
        cycle_ids = self._clusters_for_cycle_mode(mode)
        if cycle_ids.size == 0:
            return

        self.cycle_mode = mode
        self.cycle_cluster_ids = cycle_ids
        if self.selected_clust_idx not in self.cycle_cluster_ids:
            self.set_selected_cluster(int(self.cycle_cluster_ids[0]))

    def _sync_cycle_ids_with_selected(self):
        """Ensure keyboard cycling list contains the selected cluster."""
        if not hasattr(self, 'cycle_cluster_ids') or self.cycle_cluster_ids.size == 0:
            self.cycle_mode = 'all'
            self.cycle_cluster_ids = self.clust_idx

        if self.selected_clust_idx not in self.cycle_cluster_ids:
            self.cycle_mode = 'all'
            self.cycle_cluster_ids = self.clust_idx


    # --------------------------------------------------------------------------------------------
    # Cluster selection
    # --------------------------------------------------------------------------------------------

    def _update_background_for_selected_cluster(self):
        """Set central widget background to the selected cluster color."""
        rgba = self.cluster_rgba[self.selected_clust_idx].astype(int)
        self.central_widget.setStyleSheet(
            f"background-color: rgba({rgba[0]}, {rgba[1]}, {rgba[2]}, 120);"
        )

    def _update_selected_spot(self, spot):
        """Highlight only the selected spot; reset previous selection style."""
        if getattr(self, "selected_spot", None) is not None and self.selected_spot is not spot:
            self.selected_spot.setPen(DEFAULT_PEN)
            self.selected_spot.setSize(DEFAULT_SIZE)

        self.selected_spot = spot
        self.selected_spot.setPen(SELECTED_PEN)
        self.selected_spot.setSize(SELECTED_SIZE)


    def _update_plots_for_cluster(self, clust_idx: int):
        """Update all dependent plots for a given cluster index."""
        self.compute_data_for_cluster(clust_idx)
        self.plot_amplitude_distribution()
        self.plot_mean_amplitude()
        self.plot_mean_firing_rate()
        self.plot_estimated_depth()
        self.plot_amplitude_histogram()
        self.plot_firing_rate_histogram()

    def set_selected_cluster(self, clust_idx: int):
        """Select a cluster, update highlight, and refresh dependent plots."""

        self.selected_clust_idx = int(clust_idx)
        self._sync_cycle_ids_with_selected()
        clus_id = self.analyzer.sorting.get_unit_ids()[self.selected_clust_idx]
        self.fig1.setTitle(f"<span style='font-size:10pt'>Cluster {clus_id}, index {self.selected_clust_idx}</span>")
        self._update_selected_spot(self.spot_by_cluster_idx[self.selected_clust_idx])
        self._update_background_for_selected_cluster()
        self._update_plots_for_cluster(clust_idx)


def get_color(
        cmap_name: str, cbin: int = 256
    ):
        """
        Generate a pyqtgraph-compatible color map, LUT, and gradient from a given colormap.

        Parameters
        ----------
        cmap_name : str
            Name of the Matplotlib colormap.
        cbin : int, default=256
            Number of discrete bins for the LUT.

        Returns
        -------
        map : pg.ColorMap
            A pyqtgraph ColorMap object.
        lut : np.ndarray
            Lookup table for color mapping.
        grad : QtGui.QLinearGradient
            Gradient object for rendering the bar.
        """
        mpl_cmap = matplotlib.colormaps[cmap_name]
        if isinstance(mpl_cmap, mpl.colors.LinearSegmentedColormap):
            cbins = np.linspace(0.0, 1.0, cbin)
            colors = (mpl_cmap(cbins)[np.newaxis, :, :3][0]).tolist()
        else:
            colors = mpl_cmap.colors
        colors = [(np.array(c) * 255).astype(int).tolist() + [255.0] for c in colors]
        positions = np.linspace(0, 1, len(colors))
        cmap = pg.ColorMap(positions, colors)
        lut = cmap.getLookupTable()
        grad = cmap.getGradient()

        return cmap, lut, grad


if __name__ == "__main__":
    app = QApplication([])
    gui = MainGUI()
    gui.show()
    app.exec_()
