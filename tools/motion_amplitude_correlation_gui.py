"""
motion_amplitude_correlation.py

For a given session and probe, compare KS2 vs KS4 spike sorting runs by asking:
how many good/mua units have amplitudes that are statistically independent of
probe motion?

Per unit: correlate spike amplitudes with the motion estimate interpolated at
each spike's time and the unit's mean depth.

Usage
-----
python motion_amplitude_correlation.py \
    --mouse AB163 \
    --session AB163_20250419_155630 \
    --probe imec0 \
    [--base M:/analysis/Axel_Bisi/data] \
    [--alpha 0.05] \
    [--min-spikes 50]
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr

import one.alf.io as alfio
import spikeinterface.full as si
import spikeinterface.preprocessing as sp
from spikeinterface.postprocessing.unit_locations import dtype_localize_by_method
from spikeinterface.extractors.alfsortingextractor import ALFSortingExtractor
from spikeinterface.metrics.quality import ComputeQualityMetrics
from spikeinterface.sortingcomponents.tools import (
    create_sorting_analyzer_with_existing_templates,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")

KEEP_LABELS = {'GOOD', 'MUA'}
RUN_COLORS  = ['#e69f00', '#56b4e9']


# ---------------------------------------------------------------------------
# Path resolution  (mirrors SessionPickerDialog)
# ---------------------------------------------------------------------------

def resolve_runs(base, mouse, session, probe):
    session_ephys = base / mouse / session / "Ephys"
    imec_dirs = sorted(session_ephys.glob(f"catgt_*/*_{probe}"))
    if not imec_dirs:
        raise FileNotFoundError(f"No imec directory for probe {probe} in {session_ephys}")
    imec = imec_dirs[0]

    motion_path = imec / "dredge_fast" / "motion"
    if not motion_path.is_dir():
        raise FileNotFoundError(f"Motion folder not found: {motion_path}")

    runs = {}
    for ks_dir in sorted(imec.iterdir()):
        if not ks_dir.is_dir() or 'kilosort' not in ks_dir.name.lower():
            continue
        is_ks2      = ks_dir.name.lower() == "kilosort2"
        ibl_format  = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
        metrics     = (ks_dir / "cluster_bc_unitType.tsv" if is_ks2
                       else ks_dir / "sorter_output" / "cluster_bc_unitType.tsv")
        if ibl_format.is_dir() and metrics.exists():
            imec_part  = next((p for p in reversed(imec.name.split('_'))
                               if p.startswith('imec')), 'imec0')
            runs[ks_dir.name] = {
                'raw_path':    imec,
                'stream_name': f"{imec_part}.ap",
                'ibl_format':  ibl_format,
                'metrics':     metrics,
                'motion_path': motion_path,
            }

    if not runs:
        raise FileNotFoundError(f"No valid kilosort runs found under {imec}")
    return runs


# ---------------------------------------------------------------------------
# Data loading  (identical to GUI's load_data_with_spike_interface)
# ---------------------------------------------------------------------------

def load_analyzer(raw_path, stream_name, ibl_format, metrics_file):
    recording = si.read_spikeglx(folder_path=raw_path, stream_name=stream_name)
    recording.reset_times()
    fs      = recording.sampling_frequency
    sorting = ALFSortingExtractor(ibl_format, fs)

    spikes        = alfio.load_object(ibl_format, "spikes")
    templates_obj = alfio.load_object(ibl_format, "templates")
    templates_arr = np.nan_to_num(templates_obj["waveforms"])
    winds         = templates_obj["waveformsChannels"]

    sparsity = np.zeros((len(templates_arr), 384), dtype=bool)
    for i in range(len(templates_arr)):
        mask = winds[i] < 384
        sparsity[i, winds[i, mask]] = True

    templates = si.Templates(
        templates_arr * 1e6,
        sampling_frequency=fs,
        nbefore=templates_arr.shape[2] // 2,
        is_in_uV=True,
        sparsity_mask=sparsity,
        probe=recording.get_probe(),
        channel_ids=recording.channel_ids,
    )

    positions      = np.zeros(len(spikes["depths"]),
                              dtype=dtype_localize_by_method["center_of_mass"])
    positions["y"] = spikes["depths"]

    analyzer = create_sorting_analyzer_with_existing_templates(
        sorting, recording,
        spike_amplitudes=1e6 * spikes["amps"],
        spike_locations=positions,
        noise_levels=np.ones(recording.get_num_channels()),
        templates=templates,
    )

    metrics = pd.read_csv(metrics_file, sep='\t')
    metrics['gui_label'] = metrics['bc_unitType']
    analyzer.extensions["quality_metrics"] = ComputeQualityMetrics(analyzer)
    analyzer.extensions["quality_metrics"].set_data('metrics', metrics)
    analyzer.extensions["quality_metrics"].run_info["run_completed"] = True
    analyzer.extensions["quality_metrics"].run_info["runtime_s"]     = 0

    return analyzer


# ---------------------------------------------------------------------------
# Per-unit correlation
# ---------------------------------------------------------------------------

def compute_correlations(analyzer, motion, min_spikes=50):
    """
    For every unit, pull spike times / amplitudes / depths from the analyzer
    (exactly as the GUI does), interpolate motion at each spike, then run
    Pearson r against amplitude.
    """
    spike_vec  = analyzer.sorting.to_spike_vector()
    all_times  = analyzer.recording.get_times()[spike_vec["sample_index"]]
    all_amps   = analyzer.get_extension("spike_amplitudes").data['amplitudes']
    all_depths = analyzer.get_extension("spike_locations").data["spike_locations"]["y"]
    all_units  = spike_vec["unit_index"]

    qm     = analyzer.get_extension('quality_metrics').get_data()
    labels = qm['gui_label'].values   # aligned to unit_index order

    records = []
    for unit_idx in np.unique(all_units):
        mask  = all_units == unit_idx
        n     = mask.sum()
        label = labels[unit_idx] if unit_idx < len(labels) else 'noise'

        if n < min_spikes:
            continue

        t          = all_times[mask]
        a          = all_amps[mask]
        mean_depth = float(np.nanmedian(all_depths[mask]))

        # Interpolate motion at each spike's time at this unit's mean depth
        motion_at_spikes = motion.get_displacement_at_time_and_depth(
            times_s=t,
            locations_um=np.array([mean_depth]),
            segment_index=0,
            grid=False,
        ).ravel()

        valid = np.isfinite(motion_at_spikes) & np.isfinite(a)
        if valid.sum() < min_spikes:
            continue

        r, p = pearsonr(motion_at_spikes[valid], a[valid])

        records.append({
            'unit_index': int(unit_idx),
            'label':      label,
            'n_spikes':   int(n),
            'r':          float(r),
            'p_value':    float(p),
            'depth':      mean_depth,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fraction_non_significant(df, alpha):
    rows = []
    for label in ['GOOD', 'MUA', 'ALL']:
        sub = df[df['label'].isin(KEEP_LABELS)] if label == 'ALL' else df[df['label'] == label]
        rows.append({
            'label':    label,
            'fraction': (sub['p_value'] > alpha).mean() if len(sub) else np.nan,
            'n':        len(sub),
        })
    return pd.DataFrame(rows)


def plot_results(results, alpha, save_path=None):
    run_names = list(results.keys())
    colors    = RUN_COLORS[:len(run_names)]

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        f"Amplitude–motion correlation  |  α = {alpha}  |  "
        + "  vs  ".join(run_names),
        fontsize=13, fontweight='bold',
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # --- fraction of non-significant units per label ----------------------
    ax_frac = fig.add_subplot(gs[0, :2])
    fracs   = {name: fraction_non_significant(df, alpha) for name, df in results.items()}
    labels  = ['GOOD', 'MUA', 'ALL']
    x       = np.arange(len(labels))
    w       = 0.8 / len(run_names)

    for i, (name, color) in enumerate(zip(run_names, colors)):
        fdf  = fracs[name].set_index('label')
        vals = [fdf.loc[l, 'fraction'] if l in fdf.index else np.nan for l in labels]
        ns   = [fdf.loc[l, 'n']        if l in fdf.index else 0       for l in labels]
        bars = ax_frac.bar(x + i * w, vals, width=w * 0.9,
                           color=color, alpha=0.85, label=name)
        for bar, n, v in zip(bars, ns, vals):
            if np.isfinite(v):
                ax_frac.text(bar.get_x() + bar.get_width() / 2,
                             v + 0.02, f'n={n}',
                             ha='center', va='bottom', fontsize=8)

    ax_frac.set_xticks(x + w * (len(run_names) - 1) / 2)
    ax_frac.set_xticklabels(labels)
    ax_frac.set_ylim(0, 1.15)
    ax_frac.set_ylabel('Fraction non-significant  (p > α)')
    ax_frac.set_title('Non-significant units by label')
    ax_frac.legend()
    ax_frac.axhline(1.0, color='k', lw=0.5, ls='--')

    # --- GOOD vs MUA fraction scatter ------------------------------------
    ax_sc = fig.add_subplot(gs[0, 2])
    for name, color in zip(run_names, colors):
        fdf = fracs[name].set_index('label')
        fg  = fdf.loc['GOOD', 'fraction'] if 'GOOD' in fdf.index else np.nan
        fm  = fdf.loc['MUA',  'fraction'] if 'MUA'  in fdf.index else np.nan
        ax_sc.scatter(fg, fm, color=color, s=120, zorder=3)
        ax_sc.annotate(name, (fg, fm), textcoords='offset points',
                       xytext=(6, 4), fontsize=8)
    ax_sc.set_xlabel('Fraction non-sig  GOOD')
    ax_sc.set_ylabel('Fraction non-sig  MUA')
    ax_sc.set_xlim(0, 1.05); ax_sc.set_ylim(0, 1.05)
    ax_sc.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.4)
    ax_sc.set_title('GOOD vs MUA fraction')

    # --- r distributions per label ---------------------------------------
    for col_idx, label in enumerate(['GOOD', 'MUA']):
        ax = fig.add_subplot(gs[1, col_idx])
        for name, color in zip(run_names, colors):
            sub = results[name][results[name]['label'] == label]['r']
            if len(sub) == 0:
                continue
            ax.hist(sub, bins=30, range=(-1, 1), alpha=0.55,
                    color=color, label=name, density=True)
        ax.axvline(0, color='k', lw=1, ls='--')
        ax.set_xlabel('Pearson r  (amplitude ~ motion)')
        ax.set_ylabel('Density')
        ax.set_title(f'r distribution  –  {label}')
        ax.legend(fontsize=8)

    # --- log p-value distribution ----------------------------------------
    ax_p = fig.add_subplot(gs[1, 2])
    for name, color in zip(run_names, colors):
        sub = results[name][results[name]['label'].isin(KEEP_LABELS)]['p_value']
        if len(sub) == 0:
            continue
        ax_p.hist(np.log10(sub.clip(1e-10, 1)), bins=30, alpha=0.55,
                  color=color, label=name, density=True)
    ax_p.axvline(np.log10(alpha), color='k', lw=1.2, ls='--',
                 label=f'α = {alpha}')
    ax_p.set_xlabel('log₁₀(p-value)')
    ax_p.set_ylabel('Density')
    ax_p.set_title('p-value distribution  (GOOD + MUA)')
    ax_p.legend(fontsize=8)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Session picker dialog
# ---------------------------------------------------------------------------

import sys
from qtpy.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout,
    QLabel, QComboBox, QDialogButtonBox, QDoubleSpinBox, QSpinBox,
)
from qtpy.QtCore import Qt


class PickerDialog(QDialog):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Motion–Amplitude Correlation")
        self.resize(500, 300)
        self._imec_dirs = {}

        layout = QVBoxLayout(self)

        for label_text, attr in [
            ("Mouse:",   "mouse_combo"),
            ("Session:", "session_combo"),
            ("Probe:",   "probe_combo"),
        ]:
            row = QHBoxLayout()
            row.addWidget(QLabel(label_text))
            combo = QComboBox()
            setattr(self, attr, combo)
            row.addWidget(combo)
            layout.addLayout(row)

        # Alpha
        row = QHBoxLayout()
        row.addWidget(QLabel("Alpha (p-value threshold):"))
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.001, 0.5)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.05)
        self.alpha_spin.setDecimals(3)
        row.addWidget(self.alpha_spin)
        layout.addLayout(row)

        # Min spikes
        row = QHBoxLayout()
        row.addWidget(QLabel("Min spikes per unit:"))
        self.min_spikes_spin = QSpinBox()
        self.min_spikes_spin.setRange(1, 10000)
        self.min_spikes_spin.setValue(50)
        row.addWidget(self.min_spikes_spin)
        layout.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._populate_mice()
        self.mouse_combo.currentIndexChanged.connect(self._on_mouse_changed)
        self.session_combo.currentIndexChanged.connect(self._on_session_changed)
        self._on_mouse_changed()

    def _populate_mice(self):
        mice = sorted(d.name for d in BASE_DATA_PATH.iterdir()
                      if d.is_dir() and not d.name.startswith('.'))
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
        imec_dirs = sorted((BASE_DATA_PATH / mouse / session).glob("Ephys/catgt_*/*_imec*"))
        self._imec_dirs = {p.name: p for p in imec_dirs}
        self.probe_combo.addItems(list(self._imec_dirs))

    def get_params(self):
        return {
            'mouse':      self.mouse_combo.currentText(),
            'session':    self.session_combo.currentText(),
            'probe':      self.probe_combo.currentText().split('_')[-1],  # e.g. imec0
            'alpha':      self.alpha_spin.value(),
            'min_spikes': self.min_spikes_spin.value(),
        }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app = QApplication(sys.argv)

    dialog = PickerDialog()
    if dialog.exec_() != QDialog.Accepted:
        sys.exit(0)

    p = dialog.get_params()
    print(f"Resolving paths for {p['mouse']} / {p['session']} / {p['probe']} …")
    runs = resolve_runs(BASE_DATA_PATH, p['mouse'], p['session'], p['probe'])
    print(f"  runs found: {list(runs.keys())}")

    results = {}
    for run_name, paths in runs.items():
        print(f"\nLoading {run_name} …")
        analyzer = load_analyzer(
            paths['raw_path'], paths['stream_name'],
            paths['ibl_format'], paths['metrics'],
        )
        motion = sp.load_motion_info(paths['motion_path'])['motion']

        print("  Computing correlations …")
        df = compute_correlations(analyzer, motion, min_spikes=p['min_spikes'])
        results[run_name] = df

        kept = df[df['label'].isin(KEEP_LABELS)]
        frac = (kept['p_value'] > p['alpha']).mean() if len(kept) else np.nan
        print(f"  {len(df)} units  |  {len(kept)} GOOD/MUA  |  "
              f"{frac:.1%} non-significant (p > {p['alpha']})")

    plot_results(results, alpha=p['alpha'])
