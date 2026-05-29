"""
motion_amplitude_correlation_gui.py

For a given session and probe, compare KS2 vs KS4 spike sorting runs:
  - Correlate spike amplitudes with probe motion (per spike, interpolated)
  - Correlate firing rate with probe motion (2-second bins)

Both are shown on a single figure.

Data loaded directly from ALF objects — no SpikeInterface analyzer needed.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr

import one.alf.io as alfio
import spikeinterface.preprocessing as sp  # only for load_motion_info

from qtpy.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout,
    QLabel, QComboBox, QDialogButtonBox, QDoubleSpinBox, QSpinBox,
)
from qtpy.QtCore import Qt

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")

KEEP_LABELS = {'GOOD', 'MUA'}
RUN_COLORS  = ['royalblue', 'firebrick']
FR_BIN_SIZE = 2.0   # seconds


def _run_label(name):
    """Display label for a kilosort run directory name."""
    return 'ks2' if '2' in name.lower().replace('kilosort', '') else 'dredge+ks4'


# ---------------------------------------------------------------------------
# Path resolution
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
        is_ks2     = ks_dir.name.lower() == "kilosort2"
        ibl_format = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
        metrics      = (ks_dir / "cluster_bc_unitType.tsv" if is_ks2
                        else ks_dir / "sorter_output" / "cluster_bc_unitType.tsv")
        cluster_info = (ks_dir / "cluster_info.tsv" if is_ks2
                        else ks_dir / "sorter_output" / "cluster_info.tsv")
        if ibl_format.is_dir() and metrics.exists():
            runs[ks_dir.name] = {
                'ibl_format':   ibl_format,
                'metrics':      metrics,
                'cluster_info': cluster_info,
                'motion_path':  motion_path,
                'ks_dir':       ks_dir,
                'is_ks2':       is_ks2,
            }

    if not runs:
        raise FileNotFoundError(f"No valid kilosort runs found under {imec}")
    return runs


# ---------------------------------------------------------------------------
# Data loading — ALF only
# ---------------------------------------------------------------------------

def load_spikes(ibl_format, metrics_file, cluster_info_file):
    """
    Returns a dict with:
      times          : (N,) float – spike times in seconds
      amps           : (N,) float – spike amplitudes in µV
      clusters       : (N,) int   – cluster index per spike
      labels         : (C,) str   – label per cluster (indexed by cluster id)
      cluster_depths : (C,) float – peak-channel depth per cluster in µm,
                                    from cluster_info.tsv
    """
    spikes  = alfio.load_object(ibl_format, "spikes")
    metrics = pd.read_csv(metrics_file, sep='\t')
    labels  = metrics['bc_unitType'].values   # ordered by cluster_id

    ci = pd.read_csv(cluster_info_file, sep='\t')
    cluster_depths = ci['depth'].values       # one depth per cluster, ordered by cluster_id

    return {
        'times':          spikes['times'],
        'amps':           spikes['amps'] * 1e6,
        'clusters':       spikes['clusters'].astype(int),
        'labels':         labels,
        'cluster_depths': cluster_depths,
    }


# ---------------------------------------------------------------------------
# Correlation helpers
# ---------------------------------------------------------------------------

def _interp_motion(motion, times_s, depth_um):
    """
    Interpolate DREDge motion at an array of times and a single fixed depth.
    Both axes are clipped to grid bounds before querying.
    Returns an array of NaN if depth_um is NaN or non-finite.
    """
    if not np.isfinite(depth_um):
        return np.full(len(times_s), np.nan)

    t_min, t_max = motion.temporal_bins_s[0][0],  motion.temporal_bins_s[0][-1]
    d_min, d_max = motion.spatial_bins_um[0],      motion.spatial_bins_um[-1]

    times_clipped = np.clip(times_s,  t_min, t_max)
    depth_clipped = float(np.clip(depth_um, d_min, d_max))

    return motion.get_displacement_at_time_and_depth(
        times_s=times_clipped,
        locations_um=np.full(len(times_clipped), depth_clipped),
        segment_index=0,
        grid=False,
    ).ravel()


def compute_correlations(data, motion, min_spikes=50):
    """
    For every unit, Pearson r between:
      - amplitude ~ motion   (one value per spike)
      - firing rate ~ motion (2-second bins)

    Returns a DataFrame with one row per unit.
    """
    times          = data['times']
    amps           = data['amps']
    clusters       = data['clusters']
    labels         = data['labels']
    cluster_depths = data['cluster_depths']

    # FR binning over entire recording
    t_end      = times[-1]
    bin_edges  = np.arange(0.0, t_end + FR_BIN_SIZE, FR_BIN_SIZE)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    records = []
    for uid in np.unique(clusters):
        mask  = clusters == uid
        n     = mask.sum()
        label = labels[uid] if uid < len(labels) else 'noise'

        if n < min_spikes:
            continue

        t     = times[mask]
        a     = amps[mask]
        depth = (float(cluster_depths[uid]) if uid < len(cluster_depths)
                 else float('nan'))

        if not np.isfinite(depth):
            continue

        # --- amplitude ~ motion ------------------------------------------
        mot_spikes = _interp_motion(motion, t, depth)
        valid_amp  = np.isfinite(mot_spikes) & np.isfinite(a)
        if valid_amp.sum() < min_spikes:
            r_amp, p_amp = np.nan, np.nan
        else:
            r_amp, p_amp = pearsonr(mot_spikes[valid_amp], a[valid_amp])

        # --- firing rate ~ motion ----------------------------------------
        counts, _  = np.histogram(t, bins=bin_edges)
        fr         = counts / FR_BIN_SIZE             # Hz
        mot_bins   = _interp_motion(motion, bin_centers, depth)
        valid_fr   = np.isfinite(mot_bins) & np.isfinite(fr)
        if valid_fr.sum() < 10:
            r_fr, p_fr = np.nan, np.nan
        else:
            r_fr, p_fr = pearsonr(mot_bins[valid_fr], fr[valid_fr])

        records.append({
            'unit_id':   int(uid),
            'label':     label,
            'n_spikes':  int(n),
            'r_amp':     float(r_amp),
            'p_amp':     float(p_amp),
            'r_fr':      float(r_fr),
            'p_fr':      float(p_fr),
            'depth':     depth,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Presence ratio loading
# ---------------------------------------------------------------------------

def load_presence_ratio(ks_dir, is_ks2):
    """
    KS2 : <ks_dir>/qMetrics/cluster_presence_ratio.tsv
    KS4 : <ks_dir>/sorter_output/cluster_presenceRatio.tsv

    Returns a 1-D numpy array of presence ratio values, or None if not found.
    """
    path = (ks_dir / "qMetrics" / "cluster_presence_ratio.tsv" if is_ks2
            else ks_dir / "sorter_output" / "cluster_presenceRatio.tsv")
    if not path.exists():
        print(f"  [warn] presence ratio not found: {path}")
        return None
    df = pd.read_csv(path, sep='\t')
    # Identify the ratio column robustly
    col = next((c for c in df.columns if 'ratio' in c.lower() or 'presence' in c.lower()),
               df.columns[-1])
    return df[col].dropna().values


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _frac_nonsig(df, col_p, alpha):
    rows = []
    for label in ['GOOD', 'MUA', 'ALL']:
        sub = (df[df['label'].isin(KEEP_LABELS)] if label == 'ALL'
               else df[df['label'] == label])
        sub = sub.dropna(subset=[col_p])
        rows.append({
            'label':    label,
            'fraction': (sub[col_p] > alpha).mean() if len(sub) else np.nan,
            'n':        len(sub),
        })
    return pd.DataFrame(rows)


def _pick_example_units(df, col_r, col_p, alpha):
    """
    From GOOD/MUA units with valid r, return:
      uid_corr   – unit with highest |r| (most correlated)
      uid_indep  – unit with lowest |r| among those with p > alpha
    """
    sub = df[df['label'].isin(KEEP_LABELS)].dropna(subset=[col_r, col_p])
    if sub.empty:
        return None, None
    uid_corr  = sub.loc[sub[col_r].abs().idxmax(), 'unit_id']
    nonsig    = sub[sub[col_p] > alpha]
    uid_indep = nonsig.loc[nonsig[col_r].abs().idxmin(), 'unit_id'] if not nonsig.empty else None
    return int(uid_corr), (int(uid_indep) if uid_indep is not None else None)


def _scatter_data_for_unit(data, motion, uid):
    """Return (mot_spikes, amps, mot_bins, fr) for a single unit."""
    mask  = data['clusters'] == uid
    t     = data['times'][mask]
    a     = data['amps'][mask]
    depth = float(data['cluster_depths'][uid])

    mot_spikes = _interp_motion(motion, t, depth)

    t_end      = data['times'][-1]
    bin_edges  = np.arange(0.0, t_end + FR_BIN_SIZE, FR_BIN_SIZE)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts, _  = np.histogram(t, bins=bin_edges)
    fr         = counts / FR_BIN_SIZE
    mot_bins   = _interp_motion(motion, bin_centers, depth)

    return mot_spikes, a, mot_bins, fr


def _sq(ax):
    """Make an axes box square."""
    ax.set_box_aspect(1)
    return ax


def plot_results(results, run_data, run_motions, presence, alpha, save_path=None):
    """
    Single figure  (4 rows × 4 cols):
      Row 0 : fraction non-significant bar charts  (amp | FR)
      Row 1 : r-distribution histograms            (amp GOOD | amp MUA | FR GOOD | FR MUA)
      Row 2 : example scatter plots                (amp corr | amp indep | FR corr | FR indep)
      Row 3 : presence ratio comparison            (histogram spanning 2 | boxplot spanning 2)

    Bar annotations: number of units in that label group with a valid r value.
    Example units are taken from the first run that has GOOD/MUA units.
    """
    run_names  = list(results.keys())
    colors     = RUN_COLORS[:len(run_names)]
    leg_labels = [_run_label(n) for n in run_names]

    # --- pick example units from the first available run ------------------
    ex_run  = run_names[0]
    ex_df   = results[ex_run]
    ex_data = run_data[ex_run]
    ex_mot  = run_motions[ex_run]

    uid_amp_corr,  uid_amp_indep  = _pick_example_units(ex_df, 'r_amp', 'p_amp', alpha)
    uid_fr_corr,   uid_fr_indep   = _pick_example_units(ex_df, 'r_fr',  'p_fr',  alpha)

    # --- figure -----------------------------------------------------------
    fig = plt.figure(figsize=(16, 22))
    fig.suptitle(
        f"Motion correlation  |  α = {alpha}  |  "
        + "  vs  ".join(leg_labels)
        + f"  |  FR bin = {FR_BIN_SIZE} s",
        fontsize=13,
    )
    gs = gridspec.GridSpec(5, 4, figure=fig)

    # --- helpers ----------------------------------------------------------
    def _bar_panel(ax, col_p, title):
        fracs  = {n: _frac_nonsig(df, col_p, alpha) for n, df in results.items()}
        labels = ['GOOD', 'MUA', 'ALL']
        x      = np.arange(len(labels))
        w      = 0.8 / len(run_names)
        for i, (name, color, lbl) in enumerate(zip(run_names, colors, leg_labels)):
            fdf  = fracs[name].set_index('label')
            vals = [fdf.loc[l, 'fraction'] if l in fdf.index else np.nan for l in labels]
            ns   = [fdf.loc[l, 'n']        if l in fdf.index else 0       for l in labels]
            bars = ax.bar(x + i * w, vals, width=w * 0.9,
                          color=color, alpha=0.85, label=lbl)
            for bar, n, v in zip(bars, ns, vals):
                if np.isfinite(v):
                    # n = number of units with a valid r in this label group
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            v + 0.02, f'{n} units',
                            ha='center', va='bottom', fontsize=7)
        ax.set_xticks(x + w * (len(run_names) - 1) / 2)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, 1.25)
        ax.set_ylabel('Fraction non-significant  (p > α)')
        ax.set_title(title, fontsize=9)
        leg = ax.legend(fontsize=8); leg.get_frame().set_visible(False)
        ax.axhline(1.0, color='k', lw=0.5, ls='--')
        _sq(ax)

    def _hist_panel(ax, col_r, unit_label, title):
        for name, color, lbl in zip(run_names, colors, leg_labels):
            sub = results[name]
            sub = sub[sub['label'] == unit_label][col_r].dropna()
            if len(sub) == 0:
                continue
            ax.hist(sub, bins=30, range=(-1, 1), alpha=0.55,
                    color=color, label=lbl, density=True)
        ax.axvline(0, color='k', lw=1, ls='--')
        ax.set_xlabel('Pearson r', fontsize=8)
        ax.set_ylabel('Density',   fontsize=8)
        ax.set_title(title, fontsize=9)
        leg = ax.legend(fontsize=7); leg.get_frame().set_visible(False)
        _sq(ax)

    def _scatter_panel(ax, uid, col_r, col_p, x_vals, y_vals, xlabel, ylabel, title_prefix):
        if uid is None:
            ax.set_visible(False); return
        row   = ex_df[ex_df['unit_id'] == uid].iloc[0]
        r_val = row[col_r]; p_val = row[col_p]
        valid = np.isfinite(x_vals) & np.isfinite(y_vals)
        ax.scatter(x_vals[valid], y_vals[valid],
                   s=6, alpha=0.35, color='k', linewidths=0)
        ax.set_xlabel(xlabel,  fontsize=8)
        ax.set_ylabel(ylabel,  fontsize=8)
        p_str = f'{p_val:.2e}' if p_val < 0.001 else f'{p_val:.3f}'
        ax.set_title(f'{title_prefix}\nr={r_val:.2f}, p={p_str}', fontsize=8)
        _sq(ax)

    # --- Row 0: fraction bars --------------------------------------------
    _bar_panel(fig.add_subplot(gs[0, :2]), 'p_amp',
               'Amplitude ~ motion  |  fraction non-significant')
    _bar_panel(fig.add_subplot(gs[0, 2:]), 'p_fr',
               'Firing rate ~ motion  |  fraction non-significant')

    # --- Row 1: r distributions ------------------------------------------
    _hist_panel(fig.add_subplot(gs[1, 0]), 'r_amp', 'GOOD', 'r(amp ~ motion)  GOOD')
    _hist_panel(fig.add_subplot(gs[1, 1]), 'r_amp', 'MUA',  'r(amp ~ motion)  MUA')
    _hist_panel(fig.add_subplot(gs[1, 2]), 'r_fr',  'GOOD', 'r(FR ~ motion)  GOOD')
    _hist_panel(fig.add_subplot(gs[1, 3]), 'r_fr',  'MUA',  'r(FR ~ motion)  MUA')

    # --- Row 2: example scatter plots ------------------------------------
    if uid_amp_corr is not None:
        ms, a, _, _ = _scatter_data_for_unit(ex_data, ex_mot, uid_amp_corr)
        _scatter_panel(fig.add_subplot(gs[2, 0]), uid_amp_corr,
                       'r_amp', 'p_amp', ms, a,
                       'Motion (µm)', 'Amplitude (µV)',
                       f'Amp ~ motion  correlated\n(unit {uid_amp_corr})')

    if uid_amp_indep is not None:
        ms, a, _, _ = _scatter_data_for_unit(ex_data, ex_mot, uid_amp_indep)
        _scatter_panel(fig.add_subplot(gs[2, 1]), uid_amp_indep,
                       'r_amp', 'p_amp', ms, a,
                       'Motion (µm)', 'Amplitude (µV)',
                       f'Amp ~ motion  independent\n(unit {uid_amp_indep})')

    if uid_fr_corr is not None:
        _, _, mb, fr = _scatter_data_for_unit(ex_data, ex_mot, uid_fr_corr)
        _scatter_panel(fig.add_subplot(gs[2, 2]), uid_fr_corr,
                       'r_fr', 'p_fr', mb, fr,
                       'Motion (µm)', 'Firing rate (Hz)',
                       f'FR ~ motion  correlated\n(unit {uid_fr_corr})')

    if uid_fr_indep is not None:
        _, _, mb, fr = _scatter_data_for_unit(ex_data, ex_mot, uid_fr_indep)
        _scatter_panel(fig.add_subplot(gs[2, 3]), uid_fr_indep,
                       'r_fr', 'p_fr', mb, fr,
                       'Motion (µm)', 'Firing rate (Hz)',
                       f'FR ~ motion  independent\n(unit {uid_fr_indep})')

    # --- Row 3: presence ratio -------------------------------------------
    valid_runs = [n for n in run_names if presence.get(n) is not None]
    if valid_runs:
        from scipy.stats import mannwhitneyu, ks_2samp

        ax_pr_hist = _sq(fig.add_subplot(gs[3, :2]))
        ax_pr_viol = _sq(fig.add_subplot(gs[3, 2:]))

        # --- Histogram (step style, cleaner than filled bars) ------------
        for name in valid_runs:
            color = colors[run_names.index(name)]
            lbl   = leg_labels[run_names.index(name)]
            vals  = presence[name]
            ax_pr_hist.hist(vals, bins=40, range=(0, 1), density=True,
                            histtype='step', lw=1.8, color=color, label=lbl)
            ax_pr_hist.hist(vals, bins=40, range=(0, 1), density=True,
                            histtype='stepfilled', alpha=0.15, color=color)
        ax_pr_hist.set_xlabel('Presence ratio', fontsize=8)
        ax_pr_hist.set_ylabel('Density',        fontsize=8)
        ax_pr_hist.set_title('Presence ratio distribution', fontsize=9)
        ax_pr_hist.set_xlim(0, 1.05)
        ax_pr_hist.axvline(1.0, color='k', lw=0.8, ls='--', alpha=0.5)
        ax_pr_hist.spines['top'].set_visible(False)
        ax_pr_hist.spines['right'].set_visible(False)
        leg = ax_pr_hist.legend(fontsize=8, loc='upper left')
        leg.get_frame().set_visible(False)

        # --- Violin + embedded box + strip --------------------------------
        viol_data  = [presence[n] for n in valid_runs]
        viol_lbls  = [leg_labels[run_names.index(n)] for n in valid_runs]
        viol_cols  = [colors[run_names.index(n)]     for n in valid_runs]
        positions  = np.arange(1, len(valid_runs) + 1)

        parts = ax_pr_viol.violinplot(
            viol_data, positions=positions,
            showmedians=False, showextrema=False, widths=0.6,
        )
        for body, color in zip(parts['bodies'], viol_cols):
            body.set_facecolor(color)
            body.set_alpha(0.35)
            body.set_edgecolor(color)
            body.set_linewidth(1.2)

        # Narrow box overlay
        for pos, vals, color in zip(positions, viol_data, viol_cols):
            q25, med, q75 = np.percentile(vals, [25, 50, 75])
            iqr = q75 - q25
            lo  = max(vals.min(), q25 - 1.5 * iqr)
            hi  = min(vals.max(), q75 + 1.5 * iqr)
            ax_pr_viol.plot([pos, pos], [lo, hi],   color=color, lw=1.2, zorder=3)
            ax_pr_viol.plot([pos, pos], [q25, q75], color=color, lw=5,   zorder=4,
                            solid_capstyle='butt', alpha=0.7)
            ax_pr_viol.scatter(pos, med, color='white', s=30, zorder=5,
                               edgecolors=color, linewidths=1.2)

        # Jittered strip
        rng = np.random.default_rng(0)
        for pos, vals, color in zip(positions, viol_data, viol_cols):
            jitter = rng.uniform(-0.08, 0.08, size=len(vals))
            ax_pr_viol.scatter(pos + jitter, vals, s=3, alpha=0.25,
                               color=color, linewidths=0, zorder=2)

        # Significance bracket between the two runs (if exactly 2)
        if len(valid_runs) == 2:
            stat, p_mw = mannwhitneyu(viol_data[0], viol_data[1], alternative='two-sided')
            p_str = (f'p = {p_mw:.2e}' if p_mw < 0.001 else f'p = {p_mw:.3f}')
            y_br  = 1.07
            ax_pr_viol.plot([1, 1, 2, 2], [y_br - 0.02, y_br, y_br, y_br - 0.02],
                            color='k', lw=0.9)
            ax_pr_viol.text(1.5, y_br + 0.01, f'MW  {p_str}',
                            ha='center', va='bottom', fontsize=7)

        ax_pr_viol.set_xticks(positions)
        ax_pr_viol.set_xticklabels(viol_lbls, fontsize=9)
        ax_pr_viol.set_ylabel('Presence ratio', fontsize=8)
        ax_pr_viol.set_title('Presence ratio by run', fontsize=9)
        ax_pr_viol.set_ylim(-0.05, 1.18)
        ax_pr_viol.set_xlim(0.4, len(valid_runs) + 0.6)
        ax_pr_viol.spines['top'].set_visible(False)
        ax_pr_viol.spines['right'].set_visible(False)



        # --- Yield curve --------------------------------------------------
        HIGH_PR = 0.8          # threshold used for the unit-count panel
        ax_yield = _sq(fig.add_subplot(gs[4, :2]))
        thresholds = np.linspace(0, 1, 200)
        for name in valid_runs:
            color = colors[run_names.index(name)]
            lbl   = leg_labels[run_names.index(name)]
            vals  = presence[name]
            counts = np.array([(vals > thr).sum() / len(vals) for thr in thresholds])
            ax_yield.plot(thresholds, counts, color=color, lw=2, label=lbl)
        ax_yield.axvline(HIGH_PR, color='k', lw=0.8, ls='--', alpha=0.5,
                         label=f'PR > {HIGH_PR}')
        ax_yield.set_xlabel('Presence ratio threshold', fontsize=8)
        ax_yield.set_ylabel('fraction units above threshold', fontsize=8)
        ax_yield.set_title('Yield curve', fontsize=9)
        ax_yield.set_xlim(0, 1)
        ax_yield.set_ylim(bottom=0)
        ax_yield.spines['top'].set_visible(False)
        ax_yield.spines['right'].set_visible(False)
        leg = ax_yield.legend(fontsize=8); leg.get_frame().set_visible(False)

        # --- Unit count vs quality ----------------------------------------
        ax_cnt = _sq(fig.add_subplot(gs[4, 2:]))
        categories   = ['All units', f'PR > {HIGH_PR}']
        n_categories = len(categories)
        x    = np.arange(n_categories)
        w    = 0.7 / len(valid_runs)
        for i, name in enumerate(valid_runs):
            color  = colors[run_names.index(name)]
            lbl    = leg_labels[run_names.index(name)]
            vals   = presence[name]
            counts = [len(vals), int((vals > HIGH_PR).sum())]
            bars   = ax_cnt.bar(x + i * w, counts, width=w * 0.9,
                                color=color, alpha=0.85, label=lbl)
            for bar, count in zip(bars, counts):
                ax_cnt.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + ax_cnt.get_ylim()[1] * 0.01,
                            str(count), ha='center', va='bottom', fontsize=7.5)
        ax_cnt.set_xticks(x + w * (len(valid_runs) - 1) / 2)
        ax_cnt.set_xticklabels(categories, fontsize=8)
        ax_cnt.set_ylabel('# units', fontsize=8)
        ax_cnt.set_title('Total vs. high-presence units', fontsize=9)
        ax_cnt.spines['top'].set_visible(False)
        ax_cnt.spines['right'].set_visible(False)
        leg = ax_cnt.legend(fontsize=8); leg.get_frame().set_visible(False)
        # Re-annotate after ylim is known
        ax_cnt.relim(); ax_cnt.autoscale_view()

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=500, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()


# ---------------------------------------------------------------------------
# Session picker dialog
# ---------------------------------------------------------------------------

class PickerDialog(QDialog):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Motion Correlation")
        self.resize(480, 260)
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

        row = QHBoxLayout()
        row.addWidget(QLabel("Alpha (p-value threshold):"))
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.001, 0.5)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.01)
        self.alpha_spin.setDecimals(3)
        row.addWidget(self.alpha_spin)
        layout.addLayout(row)

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
        imec_dirs = sorted(
            (BASE_DATA_PATH / mouse / session).glob("Ephys/catgt_*/*_imec*")
        )
        self._imec_dirs = {p.name: p for p in imec_dirs}
        self.probe_combo.addItems(list(self._imec_dirs))

    def get_params(self):
        return {
            'mouse':      self.mouse_combo.currentText(),
            'session':    self.session_combo.currentText(),
            'probe':      self.probe_combo.currentText().split('_')[-1],
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

    results     = {}
    run_data    = {}
    run_motions = {}
    presence    = {}

    for run_name, paths in runs.items():
        print(f"\nLoading {run_name} …")
        data   = load_spikes(paths['ibl_format'], paths['metrics'], paths['cluster_info'])
        motion = sp.load_motion_info(paths['motion_path'])['motion']

        print(f"  Computing correlations ({len(np.unique(data['clusters']))} clusters) …")
        df = compute_correlations(data, motion, min_spikes=p['min_spikes'])
        results[run_name]     = df
        run_data[run_name]    = data
        run_motions[run_name] = motion

        pr = load_presence_ratio(paths['ks_dir'], paths['is_ks2'])
        presence[run_name] = pr
        if pr is not None:
            print(f"  Presence ratio: n={len(pr)}  median={np.nanmedian(pr):.3f}")

        kept = df[df['label'].isin(KEEP_LABELS)]
        print(f"  {len(df)} units pass min_spikes  |  {len(kept)} GOOD/MUA")
        for col_p, tag in [('p_amp', 'amp'), ('p_fr', 'FR')]:
            sub  = kept.dropna(subset=[col_p])
            frac = (sub[col_p] > p['alpha']).mean() if len(sub) else float('nan')
            print(f"    {tag}: {frac:.1%} non-significant (p > {p['alpha']})")

    plot_results(results, run_data, run_motions, presence, alpha=p['alpha'])
