import sys

from qtpy import QtGui
from qtpy.QtWidgets import (
    QWidget, QApplication, QSizePolicy, QMainWindow, QGridLayout,
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QDialogButtonBox, QTextEdit, QMessageBox, QSplitter,
)
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

# ---------------------------------------------------------------------------
# Paths – filled in by SessionPickerDialog at startup
# ---------------------------------------------------------------------------
BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")

RAW_PATH     = None
STREAM_NAME  = None
SORTING_PATH = None
METRICS_FILE = None
MOTION_PATH  = None

INTERVALS = [
    [0, 2000],
    [2000, 4000],
]

# ---------------------------------------------------------------------------
# pyqtgraph global config
# ---------------------------------------------------------------------------
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

GOOD_RGB    = np.array([40, 170, 70, 220],  dtype=np.uint8)
MUA_RGB     = np.array([0, 0, 255, 220],    dtype=np.uint8)
BAD_RGB     = np.array([210, 50, 50, 220],  dtype=np.uint8)
NOLABEL_RGB = np.array([70, 70, 70, 220],   dtype=np.uint8)

FONT = QtGui.QFont()
FONT.setPointSize(16)

DEFAULT_PEN   = pg.mkPen(None)
SELECTED_PEN  = pg.mkPen(255, 220, 0, width=6)
DEFAULT_SIZE  = 10
SELECTED_SIZE = 15

INTERVAL_COLORS = [
    (230, 159, 0, 70),
    (86, 180, 233, 70),
    (0, 158, 115, 70),
]

# Full-recording decimated voltage overview
VOLT_N_DISPLAY = 10000   # number of time points sampled across the full recording


# ===========================================================================
# Session picker dialog
# ===========================================================================

class SessionPickerDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Session")
        self.resize(580 * 2, 300 * 2)
        self._imec_dirs = {}

        layout = QVBoxLayout(self)

        for label_text, attr in [
            ("Mouse:",     "mouse_combo"),
            ("Session:",   "session_combo"),
            ("Probe:",     "probe_combo"),
            ("Kilosort:",  "ks_combo"),
        ]:
            row = QHBoxLayout()
            row.addWidget(QLabel(label_text))
            combo = QComboBox()
            setattr(self, attr, combo)
            row.addWidget(combo)
            layout.addLayout(row)

        self.status = QTextEdit()
        self.status.setReadOnly(True)
        self.status.setFixedHeight(280)
        layout.addWidget(self.status)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.ok_button = buttons.button(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._populate_mice()
        self.mouse_combo.currentIndexChanged.connect(self._on_mouse_changed)
        self.session_combo.currentIndexChanged.connect(self._on_session_changed)
        self.probe_combo.currentIndexChanged.connect(self._on_probe_changed)
        self.ks_combo.currentIndexChanged.connect(self._validate)
        self._on_mouse_changed()

    # ------------------------------------------------------------------
    # Population helpers
    # ------------------------------------------------------------------

    def _populate_mice(self):
        if not BASE_DATA_PATH.exists():
            self.status.setText(f"Base path not found:\n{BASE_DATA_PATH}")
            return
        mice = sorted(
            d.name for d in BASE_DATA_PATH.iterdir()
            if d.is_dir() and not d.name.startswith('.')
        )
        self.mouse_combo.addItems(mice)

    def _on_mouse_changed(self):
        self.session_combo.clear()
        mouse = self.mouse_combo.currentText()
        if not mouse:
            return
        sessions = sorted(
            (d.name for d in (BASE_DATA_PATH / mouse).iterdir() if d.is_dir()),
            reverse=True,
        )
        self.session_combo.addItems(sessions)
        self._on_session_changed()

    def _on_session_changed(self):
        self.probe_combo.clear()
        mouse, session = self.mouse_combo.currentText(), self.session_combo.currentText()
        if not mouse or not session:
            return
        imec_dirs = sorted(
            (BASE_DATA_PATH / mouse / session).glob("Ephys/catgt_*/*_imec*")
        )
        self._imec_dirs = {p.name: p for p in imec_dirs}
        self.probe_combo.addItems(list(self._imec_dirs))
        self._on_probe_changed()

    def _on_probe_changed(self):
        """Repopulate the kilosort combo with all subdirectories of the imec folder."""
        self.ks_combo.clear()
        probe_name = self.probe_combo.currentText()
        if not probe_name or probe_name not in self._imec_dirs:
            self._validate()
            return
        imec = self._imec_dirs[probe_name]
        ks_folders = sorted(d.name for d in imec.iterdir() if d.is_dir())
        ks_folders = [f for f in ks_folders if 'kilosort' in f]
        self.ks_combo.addItems(ks_folders)
        self._validate()

    # ------------------------------------------------------------------
    # Validation and path resolution
    # ------------------------------------------------------------------

    def _resolve_paths(self):
        probe_name = self.probe_combo.currentText()
        ks_name    = self.ks_combo.currentText()

        if not probe_name or probe_name not in self._imec_dirs:
            return None, ["No imec folder found under session/Ephys/catgt_*/"]
        if not ks_name:
            return None, ["No kilosort folder found inside the imec directory."]

        imec   = self._imec_dirs[probe_name]
        ks_dir = imec / ks_name

        # ibl_format lives next to the imec root for ks2, inside ks_dir for ks4+
        is_ks2          = ks_name.lower() == "kilosort2"
        ibl_format_path = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
        metrics_path    = ks_dir / "cluster_bc_unitType.tsv"
        motion_path     = imec / "dredge_fast" / "motion"

        missing = []
        if not ibl_format_path.is_dir():
            missing.append(f"  ibl_format/  →  {ibl_format_path}")
        if not metrics_path.exists():
            missing.append(f"  cluster_bc_unitType.tsv  →  {metrics_path}")
        if not motion_path.is_dir():
            missing.append(f"  dredge_fast/motion/  →  {motion_path}")
        if not list(imec.glob("*.ap.bin")):
            missing.append(f"  *.ap.bin  →  {imec}")

        if missing:
            return None, missing

        imec_part = next(
            (p for p in reversed(probe_name.split('_')) if p.startswith('imec')), 'imec0'
        )
        return {
            'raw_path':     imec,
            'stream_name':  f"{imec_part}.ap",
            'sorting_path': ibl_format_path,
            'metrics_file': metrics_path,
            'motion_path':  motion_path,
        }, []

    def _validate(self):
        paths, missing = self._resolve_paths()
        if missing:
            self.status.setStyleSheet("color: red")
            self.status.setText("Missing:\n" + "\n".join(missing))
            self.ok_button.setEnabled(False)
        else:
            self.status.setStyleSheet("color: green")
            self.status.setText(
                "All files found ✓\n"
                f"  Raw:      {paths['raw_path']}\n"
                f"  Sorting:  {paths['sorting_path']}\n"
                f"  Metrics:  {paths['metrics_file']}\n"
                f"  Motion:   {paths['motion_path']}"
            )
            self.ok_button.setEnabled(True)

    def get_paths(self):
        paths, _ = self._resolve_paths()
        return paths


# ===========================================================================
# Data loading
# ===========================================================================

def load_data_with_spike_interface():
    recording = si.read_spikeglx(folder_path=RAW_PATH, stream_name=STREAM_NAME)
    recording.reset_times()
    fs      = recording.sampling_frequency
    sorting = ALFSortingExtractor(SORTING_PATH, fs)

    spikes        = alfio.load_object(SORTING_PATH, "spikes")
    templates_obj = alfio.load_object(SORTING_PATH, "templates")
    templates_array = np.nan_to_num(templates_obj["waveforms"])
    winds = templates_obj["waveformsChannels"]

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

    positions = np.zeros(len(spikes["depths"]),
                         dtype=dtype_localize_by_method["center_of_mass"])
    positions["y"] = spikes["depths"]

    analyzer = create_sorting_analyzer_with_existing_templates(
        sorting, recording,
        spike_amplitudes=1e6 * spikes["amps"],
        spike_locations=positions,
        noise_levels=np.ones(recording.get_num_channels()),
        templates=templates,
    )

    metrics = pd.read_csv(METRICS_FILE, sep='\t')
    metrics['gui_label'] = metrics['bc_unitType']
    analyzer.extensions["quality_metrics"] = ComputeQualityMetrics(analyzer)
    analyzer.extensions["quality_metrics"].set_data('metrics', metrics)
    analyzer.extensions["quality_metrics"].run_info["run_completed"] = True
    analyzer.extensions["quality_metrics"].run_info["runtime_s"]     = 0

    motion = sp.load_motion_info(MOTION_PATH)['motion']
    return analyzer, motion


# ===========================================================================
# Main GUI
# ===========================================================================

class MainGUI(QMainWindow):

    def __init__(self):
        super().__init__()

        self.analyzer, self.motion = load_data_with_spike_interface()

        if INTERVALS is None:
            half = self.analyzer.get_total_duration() / 2
            self.intervals = [[0, half], [half, self.analyzer.get_total_duration()]]
        else:
            self.intervals = [list(iv) for iv in INTERVALS]

        self.init_gui()

        # Cluster colours
        self.compute_average_metrics_by_cluster()
        has_label = 'gui_label' in self.analyzer.get_extension('quality_metrics').get_data()
        if not has_label:
            self.cluster_rgba = np.tile(NOLABEL_RGB, (self.clust_idx.size, 1))
        else:
            self.cluster_rgba = np.tile(BAD_RGB, (self.clust_idx.size, 1))
            qm = self.analyzer.get_extension('quality_metrics').get_data()
            self.cluster_rgba[qm['gui_label'] == 'MUA']  = MUA_RGB
            self.cluster_rgba[qm['gui_label'] == 'GOOD'] = GOOD_RGB

        # Cache channel depths ordered by depth (used by voltage RMS)
        probe = self.analyzer.recording.get_probe()
        self._channel_depths = probe.contact_positions[:, 1]
        self._depth_sort     = np.argsort(self._channel_depths)
        self._sorted_depths  = self._channel_depths[self._depth_sort]

        # Compute full-recording decimated voltage (2000 single-frame reads)
        print("Computing decimated voltage overview…")
        self._compute_voltage_decimated()

        # Static plots
        self.plot_depth_amp_scatter()
        self.plot_drift_overview()
        self.plot_spike_raster()
        self.plot_drift_lines(self.fig2, self.motion.spatial_bins_um)
        self.plot_drift_lines(self.fig3, np.arange(self.motion.spatial_bins_um.size) * 20)
        self.plot_drift_intervals()

        self.spot_by_cluster_idx = {int(s.data()): s for s in self.scatter.points()}
        self.cycle_mode        = 'all'
        self.cycle_cluster_ids = self.clust_idx

        first_spot = self.scatter.points()[0]
        self.selected_clust_idx = int(first_spot.data())
        self.selected_spot      = None
        self.set_selected_cluster(self.selected_clust_idx)

    # -----------------------------------------------------------------------
    # Layout
    # -----------------------------------------------------------------------

    def init_gui(self):
        self.resize(1800, 1000)
        self.setWindowTitle('Drift QC GUI')
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        central_widget = QWidget()
        central_layout = QGridLayout(central_widget)
        self.setCentralWidget(central_widget)
        self.central_widget = central_widget

        # col 0 – depth vs amplitude scatter
        self.fig1 = pg.PlotWidget()
        self.fig1.setLabel('left', 'Depth along probe (um)')
        self.fig1.setLabel('bottom', 'Amplitude (uV)')
        self.fig1.getAxis('left').label.setFont(FONT)
        self.fig1.getAxis('bottom').label.setFont(FONT)

        # col 1 top – spike raster
        self.fig2 = pg.PlotWidget()
        self.fig2.setLabel('left', 'Depth along probe (um)')
        self.fig2.setLabel('bottom', 'Time (s)')
        self.fig2.getAxis('left').label.setFont(FONT)
        self.fig2.getAxis('bottom').label.setFont(FONT)

        # col 1 bottom – full-recording voltage RMS heatmap
        # Intentionally NOT y-linked to fig2: independent zoom/pan and free height
        self.fig_voltage = pg.PlotWidget()
        self.fig_voltage.setLabel('left', 'Depth (um)')
        self.fig_voltage.setLabel('bottom', 'Time (s)')
        self.fig_voltage.getAxis('left').label.setFont(FONT)
        self.fig_voltage.getAxis('bottom').label.setFont(FONT)
        self._volt_img = pg.ImageItem()
        self.fig_voltage.addItem(self._volt_img)
        # Dashed yellow line: selected cluster's mean depth
        self._volt_depth_line = pg.InfiniteLine(
            angle=0,
            pen=pg.mkPen(255, 220, 0, 220, width=2, style=Qt.DashLine),
        )
        self._volt_depth_line.setZValue(10)
        self.fig_voltage.addItem(self._volt_depth_line)

        # QSplitter in col 1: user can drag the divider to resize both panels freely
        col1_splitter = QSplitter(Qt.Vertical)
        col1_splitter.addWidget(self.fig2)
        col1_splitter.addWidget(self.fig_voltage)
        col1_splitter.setStretchFactor(0, 3)   # initial ratio: raster ~60 %
        col1_splitter.setStretchFactor(1, 2)   # voltage RMS ~40 %

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
        self.fig4.setLabel('left', 'Amplitude (uV)')
        self.fig4.getAxis('left').label.setFont(FONT)

        self.fig5 = pg.PlotItem()
        self.fr_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2))
        self.fig5.addItem(self.fr_curve)
        self.fig5.getAxis('left').setWidth(max_width)
        self.fig5.setTitle("<span style='font-size:8pt'>Mean cluster firing rate across session</span>")
        self.fig5.setLabel('left', 'Mean firing rate (spikes/s)')
        self.fig5.getAxis('left').label.setFont(FONT)

        self.fig6 = pg.PlotItem()
        self.amp_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2))
        self.fig6.addItem(self.amp_curve)
        self.fig6.getAxis('left').setWidth(max_width)
        self.fig6.setTitle("<span style='font-size:8pt'>Mean amplitude across session</span>")
        self.fig6.setXRange(0, self.motion.temporal_bins_s[0].max() + 100)
        self.fig6.setLabel('left', 'Mean amplitude (uV)')
        self.fig6.getAxis('left').label.setFont(FONT)

        self.fig7 = pg.PlotItem()
        self.drift_curve = pg.PlotCurveItem(pen=pg.mkPen(0, 0, 0, width=2))
        self.fig7.addItem(self.drift_curve)
        self.fig7.getAxis('left').setWidth(max_width)
        self.fig7.setTitle(
            "<span style='font-size:8pt'>Estimated drift at cluster depth "
            "– drag coloured regions to change intervals</span>"
        )
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

        layout      = pg.GraphicsLayout()
        widget      = pg.GraphicsLayoutWidget()
        hist_layout = pg.GraphicsLayout()
        hist_widget = pg.GraphicsLayoutWidget()

        layout.addItem(self.fig3, row=0, col=0, rowspan=3)
        layout.addItem(self.fig4, row=4, col=0, rowspan=3)
        layout.addItem(self.fig5, row=7, col=0)
        layout.addItem(self.fig6, row=8, col=0)
        layout.addItem(self.fig7, row=9, col=0)

        # Link all x-axes in col 2 to fig7 so time ranges stay identical.
        # Suppress bottom tick labels on intermediate plots so only fig7
        # (the bottom panel) labels the time axis — this also removes the
        # extra height those labels would add, keeping plot areas pixel-aligned.
        for fig in [self.fig3, self.fig4, self.fig5, self.fig6]:
            fig.setXLink(self.fig7)
            fig.getAxis('bottom').setStyle(showValues=False)

        hist_layout.addItem(self.fig8, row=0, col=0)
        hist_layout.addItem(self.fig9, row=1, col=0)

        widget.addItem(layout)
        hist_widget.addItem(hist_layout)

        central_layout.addWidget(self.fig1,       0, 0, 5, 1)
        central_layout.addWidget(col1_splitter,   0, 1, 5, 1)
        central_layout.addWidget(widget,          0, 2, 5, 1)
        central_layout.addWidget(hist_widget,     0, 3, 5, 1)

        # col 1 (raster + voltage) gets 3× more width than col 0 (scatter)
        central_layout.setColumnStretch(0, 1)
        central_layout.setColumnStretch(1, 3)
        central_layout.setColumnStretch(2, 2)
        central_layout.setColumnStretch(3, 1)

    # -----------------------------------------------------------------------
    # Data
    # -----------------------------------------------------------------------

    def compute_average_metrics_by_cluster(self):
        df = pd.DataFrame()
        df['clusters'] = self.analyzer.sorting.to_spike_vector()["unit_index"]
        df['amps']     = self.analyzer.get_extension("spike_amplitudes").data['amplitudes']
        df['depths']   = self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"]
        avgs = df.groupby('clusters').agg(['mean', 'count'])
        self.clust_idx = avgs.index.values
        self.avg_amp   = avgs['amps']['mean'].values
        self.avg_fr    = avgs['depths']['count'].values / self.analyzer.get_total_duration()
        self.avg_depth = avgs['depths']['mean'].values
        del df, avgs

    def compute_data_for_cluster(self, clust_idx: int):
        idx = self.analyzer.sorting.to_spike_vector()["unit_index"] == clust_idx
        self.clust_times = self.analyzer.recording.get_times()[
            self.analyzer.sorting.to_spike_vector()["sample_index"][idx]
        ]
        self.clust_amps = self.analyzer.get_extension("spike_amplitudes").data['amplitudes'][idx]

        n_bins      = self.motion.temporal_bins_s[0].size
        bin_indices = np.clip(
            np.searchsorted(self.motion.temporal_bin_edges_s[0], self.clust_times, side="right") - 1,
            0, n_bins - 1,
        )
        amp_sum   = np.bincount(bin_indices, weights=self.clust_amps, minlength=n_bins)
        amp_count = np.bincount(bin_indices, minlength=n_bins)

        self.mean_amp = np.full(n_bins, np.nan)
        self.mean_fr  = np.full(n_bins, np.nan)
        valid = amp_count > 0
        self.mean_amp[valid] = amp_sum[valid] / amp_count[valid]
        self.mean_fr[valid]  = amp_count[valid]

    # -----------------------------------------------------------------------
    # Static plots
    # -----------------------------------------------------------------------

    def plot_depth_amp_scatter(self):
        self.scatter = pg.ScatterPlotItem()
        self.scatter.sigClicked.connect(self._on_plot_click)
        self.fig1.addItem(self.scatter)
        self.scatter.setData(
            x=self.avg_amp, y=self.avg_depth, data=self.clust_idx,
            brush=[pg.mkBrush(*rgba) for rgba in self.cluster_rgba],
            pen=DEFAULT_PEN, size=DEFAULT_SIZE, pxMode=True,
        )
        self.fig1.setXRange(-50, np.nanmin([np.nanmax(self.avg_amp) + 200, 1200]))
        self.fig1.setYRange(0, 4000)

    def plot_spike_raster(self):
        isnan  = ~np.isnan(
            self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"]
        )
        depths = self.analyzer.get_extension("spike_locations").data["spike_locations"]["y"][isnan]
        times  = self.analyzer.recording.get_times()[
            self.analyzer.sorting.to_spike_vector()["sample_index"][isnan]
        ]
        fr, times, depths = bincount2D(
            times, depths, xbin=0.05, ybin=10,
            ylim=[min(0, depths.min()), max(3840, depths.max())],
        )
        fr = fr.T
        xscale = (times[-1] - times[0]) / fr.shape[0]
        yscale = (depths[-1] - depths[0]) / fr.shape[1]
        levels = np.quantile(np.mean(fr, axis=0), [0, 1])
        fr_img = pg.ImageItem()
        fr_img.setImage(fr)
        fr_img.setTransform(QTransform(xscale, 0, 0, 0, yscale, 0, times[0], depths[0], 1))
        fr_img.setOpacity(0.8)
        _, lut, _ = get_color('binary')
        fr_img.setLookupTable(lut)
        fr_img.setLevels(levels)
        self.fig2.addItem(fr_img)
        self.fig2.setXRange(0, self.analyzer.get_total_duration() + 50)
        self.fig2.setYRange(0, 4000)

    def plot_drift_overview(self):
        drift_interp = self.motion.get_displacement_at_time_and_depth(
            times_s=self.motion.temporal_bins_s[0],
            locations_um=np.arange(0, 3840, 40),
            segment_index=0, grid=True,
        ).T
        xscale  = (self.motion.temporal_bin_edges_s[0].max() -
                   self.motion.temporal_bin_edges_s[0].min()) / self.motion.temporal_bins_s[0].size
        yscale  = 3840 / drift_interp.shape[1]
        xoffset = self.motion.temporal_bin_edges_s[0].min()
        drift_img = pg.ImageItem()
        drift_img.setImage(drift_interp)
        drift_img.setTransform(QTransform(xscale, 0, 0, 0, yscale, 0, xoffset, 0, 1))
        _, lut, _ = get_color('seismic')
        drift_img.setLookupTable(lut)
        drift_img.setLevels((-30, 30))
        self.fig2.addItem(drift_img)

    def plot_drift_lines(self, fig, offsets):
        for disp, offset in zip(self.motion.displacement[0].T, offsets):
            curve = pg.PlotCurveItem()
            curve.setData(x=self.motion.temporal_bins_s[0], y=disp + offset, pen='r', linewidth=12)
            fig.addItem(curve)
        fig.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_drift_intervals(self):
        self._interval_regions = []
        t_min = float(np.nanmin(self.motion.temporal_bins_s[0]))
        t_max = float(np.nanmax(self.motion.temporal_bins_s[0]))
        for (start, end), color in zip(self.intervals, INTERVAL_COLORS):
            lo = max(float(start), t_min)
            hi = min(float(end),   t_max)
            if hi <= lo:
                continue
            region = pg.LinearRegionItem(
                values=(lo, hi),
                orientation='vertical',
                movable=True,
                brush=pg.mkBrush(*color),
                pen=pg.mkPen(color[0], color[1], color[2], 200),
                swapMode='sort',
            )
            region.setZValue(-20)
            region.sigRegionChanged.connect(self._on_interval_changed)
            self.fig7.addItem(region)
            self._interval_regions.append(region)

    def _on_interval_changed(self):
        self.intervals = [list(r.getRegion()) for r in self._interval_regions]
        if hasattr(self, 'mean_amp') and self.mean_amp is not None:
            self.plot_amplitude_histogram()
            self.plot_firing_rate_histogram()

    # -----------------------------------------------------------------------
    # Full-recording decimated voltage overview
    # -----------------------------------------------------------------------

    def _compute_voltage_decimated(self):
        """Stride-subsample the raw recording to VOLT_N_DISPLAY time points.
        Stores the full (n_time, n_channels_depth_sorted) array; slicing to
        50 channels is done cheaply per cluster in _update_voltage_panel."""
        rec     = self.analyzer.recording
        n_total = rec.get_num_samples()
        n_ch    = rec.get_num_channels()
        stride  = max(1, n_total // VOLT_N_DISPLAY)

        volt = np.zeros((VOLT_N_DISPLAY, n_ch), dtype=np.float32)
        for i in range(VOLT_N_DISPLAY):
            frame = i * stride
            if frame >= n_total:
                break
            volt[i] = rec.get_traces(
                start_frame=frame, end_frame=frame + 1, return_scaled=True
            )[0]

        self._voltage_dec        = volt[:, self._depth_sort]   # depth-sorted columns
        self._voltage_dec_times  = np.arange(VOLT_N_DISPLAY) * stride / rec.sampling_frequency
        self._voltage_dec_stride = stride

    def _update_voltage_panel(self):
        """Slice the cached array to 50 channels around the selected cluster
        and repaint the heatmap.  No disk I/O — pure numpy indexing."""
        N_SHOW = 50
        target = self.avg_depth[self.selected_clust_idx]

        # 50 nearest channels in depth-sorted order
        sel = np.sort(np.argsort(np.abs(self._sorted_depths - target))[:N_SHOW])
        ch_depths  = self._sorted_depths[sel]
        volt_slice = self._voltage_dec[:, sel]   # (n_time, N_SHOW)

        t0 = float(self._voltage_dec_times[0])
        t1 = float(self._voltage_dec_times[-1])
        d0 = float(ch_depths[0])
        d1 = float(ch_depths[-1])
        xscale = (t1 - t0) / max(volt_slice.shape[0] - 1, 1)
        yscale = (d1 - d0) / max(volt_slice.shape[1] - 1, 1)

        clim = float(np.nanpercentile(np.abs(volt_slice), 99))
        clim = max(clim, 1.0)

        self._volt_img.setImage(volt_slice)
        self._volt_img.setTransform(QTransform(xscale, 0, 0, 0, yscale, 0, t0, d0, 1))
        _, lut, _ = get_color('seismic')
        self._volt_img.setLookupTable(lut)
        self._volt_img.setLevels((-clim, clim))

        self._volt_depth_line.setPos(target)

        spacing = float(np.median(np.diff(ch_depths))) if len(ch_depths) > 1 else 20.0
        pad = max(spacing, 20.0) * 0.6
        self.fig_voltage.setYRange(d0 - pad, d1 + pad)
        self.fig_voltage.setXRange(0, self.analyzer.get_total_duration() + 50)

        clus_id   = self.analyzer.sorting.get_unit_ids()[self.selected_clust_idx]
        stride_s  = self._voltage_dec_stride / self.analyzer.recording.sampling_frequency
        self.fig_voltage.setTitle(
            f"<span style='font-size:8pt'>Raw voltage – cluster {clus_id} "
            f"| {N_SHOW} ch @ {target:.0f} \u00b5m "
            f"| 1 sample / {stride_s:.1f} s "
            f"| \u00b1{clim:.1f} \u00b5V</span>"
        )

    # -----------------------------------------------------------------------
    # Per-cluster plots
    # -----------------------------------------------------------------------

    def plot_amplitude_distribution(self):
        self.amp_img.clear()
        bins = [400, 200]
        H, xedges, yedges = np.histogram2d(self.clust_times, self.clust_amps * 1e6, bins=bins)
        H[H == 0] = np.nan
        xscale = (xedges[-1] - xedges[0]) / bins[0]
        yscale = (yedges[-1] - yedges[0]) / bins[1]
        self.amp_img.setImage(H)
        self.amp_img.setTransform(
            QTransform(xscale, 0, 0, 0, yscale, 0, xedges[0], yedges[0], 1)
        )
        self.amp_img.setLookupTable(pg.colormap.get("viridis").getLookupTable())
        self.fig4.setYRange(0, yedges[-1] + 20)
        self.fig4.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_mean_amplitude(self):
        self.amp_curve.setData(self.motion.temporal_bins_s[0], self.mean_amp * 1e6)
        self.fig6.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_mean_firing_rate(self):
        self.fr_curve.setData(self.motion.temporal_bins_s[0], self.mean_fr)
        self.fig5.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_estimated_depth(self):
        drift_at_depth = self.motion.get_displacement_at_time_and_depth(
            times_s=self.motion.temporal_bins_s[0],
            locations_um=[self.avg_depth[self.selected_clust_idx]],
            grid=True,
        )[0]
        self.drift_curve.setData(self.motion.temporal_bins_s[0], drift_at_depth)
        self.fig7.setXRange(0, self.analyzer.get_total_duration() + 50)

    def plot_histogram(self, data, fig, hist_items):
        for item in hist_items:
            fig.removeItem(item)
        hist_items = []

        data_max = float(np.nanmax(data)) if np.any(np.isfinite(data)) else 1.0
        if data_max <= 0:
            data_max = 1.0
        bin_edges = np.linspace(0.0, data_max * 1.05, 24)
        bin_width = np.diff(bin_edges)
        centers   = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        times     = self.motion.temporal_bins_s[0]

        for (start, end), color in zip(self.intervals, INTERVAL_COLORS):
            vals = data[(times >= start) & (times < end)]
            if vals.size == 0:
                continue
            counts, _ = np.histogram(vals, bins=bin_edges)
            total = counts.sum()
            if total == 0:
                continue
            counts = counts.astype(float) / (total * bin_width)
            bar = pg.BarGraphItem(
                x=centers, height=counts, width=bin_width * 0.92,
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

    # -----------------------------------------------------------------------
    # Interactions
    # -----------------------------------------------------------------------

    def _on_plot_click(self, _, point):
        spots = self.scatter.pointsAt(point[0].pos())
        if spots:
            self.set_selected_cluster(int(spots[0].data()))

    def keyPressEvent(self, event):
        if event.modifiers() & Qt.ShiftModifier:
            mode_map = {Qt.Key_G: 'good', Qt.Key_M: 'mua',
                        Qt.Key_B: 'bad',  Qt.Key_A: 'all'}
            mode = mode_map.get(event.key())
            if mode is not None:
                self._set_cycle_mode(mode)
                event.accept()
                return

        if event.key() in (Qt.Key_Left, Qt.Key_Right):
            self._sync_cycle_ids_with_selected()
            if self.cycle_cluster_ids.size == 0:
                event.accept()
                return
            pos  = int(np.where(self.cycle_cluster_ids == self.selected_clust_idx)[0][0])
            step = -1 if event.key() == Qt.Key_Left else 1
            self.set_selected_cluster(
                int(self.cycle_cluster_ids[(pos + step) % self.cycle_cluster_ids.size])
            )
            event.accept()
            return

        super().keyPressEvent(event)

    def _clusters_for_cycle_mode(self, mode):
        if mode == 'good': return self.clust_idx[np.all(self.cluster_rgba == GOOD_RGB, axis=1)]
        if mode == 'mua':  return self.clust_idx[np.all(self.cluster_rgba == MUA_RGB,  axis=1)]
        if mode == 'bad':  return self.clust_idx[np.all(self.cluster_rgba == BAD_RGB,  axis=1)]
        return self.clust_idx

    def _set_cycle_mode(self, mode):
        ids = self._clusters_for_cycle_mode(mode)
        if ids.size == 0:
            return
        self.cycle_mode = mode
        self.cycle_cluster_ids = ids
        if self.selected_clust_idx not in ids:
            self.set_selected_cluster(int(ids[0]))

    def _sync_cycle_ids_with_selected(self):
        if not hasattr(self, 'cycle_cluster_ids') or self.cycle_cluster_ids.size == 0:
            self.cycle_mode = 'all'
            self.cycle_cluster_ids = self.clust_idx
        if self.selected_clust_idx not in self.cycle_cluster_ids:
            self.cycle_mode = 'all'
            self.cycle_cluster_ids = self.clust_idx

    def _update_background_for_selected_cluster(self):
        rgba = self.cluster_rgba[self.selected_clust_idx].astype(int)
        self.central_widget.setStyleSheet(
            f"background-color: rgba({rgba[0]},{rgba[1]},{rgba[2]},120);"
        )

    def _update_selected_spot(self, spot):
        if getattr(self, 'selected_spot', None) is not None and self.selected_spot is not spot:
            self.selected_spot.setPen(DEFAULT_PEN)
            self.selected_spot.setSize(DEFAULT_SIZE)
        self.selected_spot = spot
        self.selected_spot.setPen(SELECTED_PEN)
        self.selected_spot.setSize(SELECTED_SIZE)

    def _update_plots_for_cluster(self, clust_idx):
        self.compute_data_for_cluster(clust_idx)
        self.plot_amplitude_distribution()
        self.plot_mean_amplitude()
        self.plot_mean_firing_rate()
        self.plot_estimated_depth()
        self.plot_amplitude_histogram()
        self.plot_firing_rate_histogram()
        self._update_voltage_panel()

    def set_selected_cluster(self, clust_idx: int):
        self.selected_clust_idx = int(clust_idx)
        self._sync_cycle_ids_with_selected()
        clus_id = self.analyzer.sorting.get_unit_ids()[self.selected_clust_idx]
        self.fig1.setTitle(
            f"<span style='font-size:10pt'>Cluster {clus_id}, "
            f"index {self.selected_clust_idx}</span>"
        )
        self._update_selected_spot(self.spot_by_cluster_idx[self.selected_clust_idx])
        self._update_background_for_selected_cluster()
        self._update_plots_for_cluster(clust_idx)


# ===========================================================================
# Utilities
# ===========================================================================

def get_color(cmap_name: str, cbin: int = 256):
    mpl_cmap = matplotlib.colormaps[cmap_name]
    if isinstance(mpl_cmap, mpl.colors.LinearSegmentedColormap):
        cbins  = np.linspace(0.0, 1.0, cbin)
        colors = mpl_cmap(cbins)[:, :3].tolist()
    else:
        colors = mpl_cmap.colors
    colors    = [(np.array(c) * 255).astype(int).tolist() + [255] for c in colors]
    positions = np.linspace(0, 1, len(colors))
    cmap      = pg.ColorMap(positions, colors)
    return cmap, cmap.getLookupTable(), cmap.getGradient()


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    app = QApplication(sys.argv)
    font = app.font()
    font.setPointSize(16)
    app.setFont(font)

    dialog = SessionPickerDialog()
    if dialog.exec_() != QDialog.Accepted:
        sys.exit(0)

    paths = dialog.get_paths()
    if paths is None:
        QMessageBox.critical(None, "Error", "Could not resolve session paths.")
        sys.exit(1)

    RAW_PATH     = paths['raw_path']
    STREAM_NAME  = paths['stream_name']
    SORTING_PATH = paths['sorting_path']
    METRICS_FILE = paths['metrics_file']
    MOTION_PATH  = paths['motion_path']

    gui = MainGUI()
    gui.show()
    sys.exit(app.exec_())
