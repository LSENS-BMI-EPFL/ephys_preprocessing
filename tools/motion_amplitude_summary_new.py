"""
motion_amplitude_summary.py

Dataset-wide amplitude–motion correlation summary.
Optimised: no raw recording loaded (ALF only), sessions with both KS2+KS4
only, parallel session processing.

Usage
-----
python motion_amplitude_summary.py [--base M:/analysis/Axel_Bisi/data]
                                   [--alpha 0.05] [--min-spikes 50]
                                   [--jobs 4]
                                   [--save summary.png] [--csv summary.csv]
"""

import argparse
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
import one.alf.io as alfio
import spikeinterface.preprocessing as sp

BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")
KEEP_LABELS    = {'GOOD', 'MUA'}
RUN_COLORS     = ['#e69f00', '#56b4e9']


# ---------------------------------------------------------------------------
# Lightweight data loading  (no recording, no SpikeInterface extractor)
# ---------------------------------------------------------------------------

def load_alf_spikes(ibl_format, metrics_file):
    """
    Load spikes and per-unit labels straight from ALF files.
    Returns arrays aligned to spikes: times, amps, depths, clusters,
    and a dict {cluster_id: label}.
    """
    spikes   = alfio.load_object(ibl_format, "spikes")
    times    = spikes['times']
    amps     = spikes['amps'] * 1e6     # → µV
    depths   = spikes['depths']
    clusters = spikes['clusters']

    metrics  = pd.read_csv(metrics_file, sep='\t')
    id_col   = next((c for c in metrics.columns
                     if c.lower() in ('cluster_id', 'cluster_ids', 'id')), None)
    if id_col is not None:
        label_map = dict(zip(metrics[id_col], metrics['bc_unitType']))
    else:
        label_map = dict(enumerate(metrics['bc_unitType']))

    return times, amps, depths, clusters, label_map


# ---------------------------------------------------------------------------
# Per-unit correlation
# ---------------------------------------------------------------------------

def compute_correlations_alf(ibl_format, metrics_file, motion, min_spikes):
    times, amps, depths, clusters, label_map = load_alf_spikes(ibl_format, metrics_file)

    records = []
    for uid in np.unique(clusters):
        mask  = clusters == uid
        if mask.sum() < min_spikes:
            continue
        label = label_map.get(int(uid), 'noise')

        t          = times[mask]
        a          = amps[mask]
        mean_depth = float(np.nanmedian(depths[mask]))

        motion_vals = motion.get_displacement_at_time_and_depth(
            times_s=t,
            locations_um=np.full(t.shape, mean_depth),
            segment_index=0,
            grid=False,
        ).ravel()

        valid = np.isfinite(motion_vals) & np.isfinite(a)
        if valid.sum() < min_spikes:
            continue

        r, p = pearsonr(motion_vals[valid], a[valid])
        records.append({'unit_id': int(uid), 'label': label,
                        'n_spikes': int(mask.sum()), 'r': float(r), 'p_value': float(p)})

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Session/probe discovery — KS2 + KS4 pairs only
# ---------------------------------------------------------------------------

KS2_NAME = 'kilosort2'

def find_valid_probes(base):
    """
    Return list of probe dicts that have BOTH a valid KS2 and KS4 run.
    """
    probes = []
    for mouse_dir in sorted(base.iterdir()):
        if not mouse_dir.is_dir() or mouse_dir.name.startswith('.'):
            continue
        for session_dir in sorted(mouse_dir.iterdir()):
            if not session_dir.is_dir():
                continue
            ephys = session_dir / "Ephys"
            if not ephys.exists():
                continue
            for imec in sorted(ephys.glob("catgt_*/*_imec*")):
                if not imec.is_dir():
                    continue
                motion_path = imec / "dredge_fast" / "motion"
                if not motion_path.is_dir():
                    continue

                ks2, ks4 = None, None
                for ks_dir in sorted(imec.iterdir()):
                    if not ks_dir.is_dir() or 'kilosort' not in ks_dir.name.lower():
                        continue
                    is_ks2 = ks_dir.name.lower() == KS2_NAME
                    ibl    = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
                    mfile  = (ks_dir / "cluster_bc_unitType.tsv" if is_ks2
                              else ks_dir / "sorter_output" / "cluster_bc_unitType.tsv")
                    if not (ibl.is_dir() and mfile.exists()):
                        continue
                    entry = {'name': ks_dir.name, 'ibl_format': ibl, 'metrics': mfile}
                    if is_ks2:
                        ks2 = entry
                    elif 'kilosort' in ks_dir.name.lower():
                        ks4 = entry   # last valid non-KS2 kilosort wins

                if ks2 and ks4:
                    probes.append({
                        'mouse':       mouse_dir.name,
                        'session':     session_dir.name,
                        'probe':       imec.name,
                        'motion_path': motion_path,
                        'ks2':         ks2,
                        'ks4':         ks4,
                    })

    return probes


# ---------------------------------------------------------------------------
# Worker (runs in subprocess)
# ---------------------------------------------------------------------------

def _process_probe(probe, alpha, min_spikes):
    tag = f"{probe['mouse']} / {probe['session']} / {probe['probe']}"
    records = []
    try:
        motion = sp.load_motion_info(probe['motion_path'])['motion']
        for key in ('ks2', 'ks4'):
            p      = probe[key]
            df     = compute_correlations_alf(p['ibl_format'], p['metrics'],
                                              motion, min_spikes)
            for label in ['GOOD', 'MUA']:
                sub = df[df['label'] == label]
                n   = len(sub)
                ns  = int((sub['p_value'] > alpha).sum()) if n else 0
                records.append({
                    'mouse':           probe['mouse'],
                    'session':         probe['session'],
                    'probe':           probe['probe'],
                    'run':             p['name'],
                    'label':           label,
                    'n_units':         n,
                    'n_nonsig':        ns,
                    'fraction_nonsig': ns / n if n else np.nan,
                })
        print(f"  ✓  {tag}")
    except Exception:
        print(f"  ✗  {tag}")
        traceback.print_exc()
    return records


# ---------------------------------------------------------------------------
# Dataset runner
# ---------------------------------------------------------------------------

def run_dataset(base, alpha, min_spikes, n_jobs):
    probes = find_valid_probes(base)
    print(f"Found {len(probes)} probes with both KS2 + KS4\n")

    all_records = []
    if n_jobs == 1:
        for probe in probes:
            all_records.extend(_process_probe(probe, alpha, min_spikes))
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = {ex.submit(_process_probe, p, alpha, min_spikes): p
                       for p in probes}
            for fut in as_completed(futures):
                all_records.extend(fut.result())

    return pd.DataFrame(all_records)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_summary(results, alpha, save_path=None):
    run_names = sorted(results['run'].unique())
    colors    = {r: RUN_COLORS[i % len(RUN_COLORS)] for i, r in enumerate(run_names)}

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(
        f"Amplitude–motion correlation  |  dataset summary  |  α = {alpha}",
        fontsize=13, fontweight='bold',
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    for col, label in enumerate(['GOOD', 'MUA']):
        ax  = fig.add_subplot(gs[0, col])
        sub = results[results['label'] == label]
        for _, row in sub[['mouse', 'session', 'probe']].drop_duplicates().iterrows():
            sd = sub[(sub['mouse'] == row['mouse']) &
                     (sub['session'] == row['session']) &
                     (sub['probe'] == row['probe'])].set_index('run')
            if all(r in sd.index for r in run_names[:2]):
                ax.plot([0, 1],
                        [sd.loc[run_names[0], 'fraction_nonsig'],
                         sd.loc[run_names[1], 'fraction_nonsig']],
                        color='#888', lw=0.8, alpha=0.4, zorder=1)
        for i, run in enumerate(run_names):
            vals   = sub[sub['run'] == run]['fraction_nonsig'].dropna()
            jitter = np.random.uniform(-0.04, 0.04, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals,
                       color=colors[run], s=28, alpha=0.75, zorder=2, label=run)
            ax.plot([i - 0.15, i + 0.15], [vals.median(), vals.median()],
                    color=colors[run], lw=2.5, zorder=3)
        ax.set_xticks(range(len(run_names)))
        ax.set_xticklabels(run_names, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel('Fraction non-significant')
        ax.set_title(f'{label} units  (median bar)')
        ax.axhline(0.5, color='k', lw=0.5, ls='--', alpha=0.3)
        if col == 0:
            ax.legend(fontsize=8)

    ax_diff = fig.add_subplot(gs[0, 2])
    if len(run_names) == 2:
        for label, marker, color in [('GOOD', 'o', '#28aa46'), ('MUA', 's', '#0000dc')]:
            sub   = results[results['label'] == label]
            diffs = []
            for _, row in sub[['mouse', 'session', 'probe']].drop_duplicates().iterrows():
                sd = sub[(sub['mouse'] == row['mouse']) &
                         (sub['session'] == row['session']) &
                         (sub['probe'] == row['probe'])].set_index('run')
                if all(r in sd.index for r in run_names[:2]):
                    diffs.append(sd.loc[run_names[1], 'fraction_nonsig'] -
                                 sd.loc[run_names[0], 'fraction_nonsig'])
            if diffs:
                jitter = np.random.uniform(-0.08, 0.08, size=len(diffs))
                ax_diff.scatter(jitter, diffs, color=color, marker=marker,
                                s=30, alpha=0.75, label=label)
                ax_diff.plot([-0.15, 0.15], [np.median(diffs)] * 2,
                             color=color, lw=2.5)
        ax_diff.axhline(0, color='k', lw=1, ls='--')
        ax_diff.set_xticks([])
        ax_diff.set_ylabel(f'Δ fraction  ({run_names[1]} − {run_names[0]})')
        ax_diff.set_title('Paired difference per session')
        ax_diff.legend(fontsize=8)
        ax_diff.set_xlim(-0.4, 0.4)

    ax_bar  = fig.add_subplot(gs[1, :2])
    labels  = ['GOOD', 'MUA']
    x, w    = np.arange(len(labels)), 0.35
    for i, run in enumerate(run_names):
        means = [results[(results['run'] == run) &
                         (results['label'] == l)]['fraction_nonsig'].dropna().mean()
                 for l in labels]
        sems  = [results[(results['run'] == run) &
                         (results['label'] == l)]['fraction_nonsig'].dropna().sem()
                 for l in labels]
        ax_bar.bar(x + i * w, means, width=w * 0.9, yerr=sems,
                   color=colors[run], alpha=0.85, label=run,
                   error_kw=dict(lw=1, capsize=4))
    ax_bar.set_xticks(x + w * (len(run_names) - 1) / 2)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylim(0, 1.15)
    ax_bar.set_ylabel('Mean fraction non-significant  (± SEM)')
    ax_bar.set_title('Dataset mean  (error bars = SEM across sessions)')
    ax_bar.legend()
    ax_bar.axhline(1.0, color='k', lw=0.5, ls='--')

    ax_tbl = fig.add_subplot(gs[1, 2])
    ax_tbl.axis('off')
    rows = [[run, label,
             f"{len(results[(results['run']==run)&(results['label']==label)])}",
             f"{results[(results['run']==run)&(results['label']==label)]['n_units'].sum():.0f}",
             f"{results[(results['run']==run)&(results['label']==label)]['fraction_nonsig'].mean():.2f}"]
            for run in run_names for label in ['GOOD', 'MUA']]
    tbl = ax_tbl.table(cellText=rows,
                       colLabels=['Run', 'Label', 'Probes', 'Units', 'Mean frac'],
                       loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.4)
    ax_tbl.set_title('Summary table', pad=12)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nFigure saved to {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--base',       default=str(BASE_DATA_PATH))
    parser.add_argument('--alpha',      type=float, default=0.05)
    parser.add_argument('--min-spikes', type=int,   default=50)
    parser.add_argument('--jobs',       type=int,   default=1,
                        help='Parallel workers (default 1; use 4-8 on a server)')
    parser.add_argument('--save',       default=None)
    parser.add_argument('--csv',        default=None)
    args = parser.parse_args()

    base = Path(args.base)
    print(f"Scanning {base} …")

    results = run_dataset(base, alpha=args.alpha,
                          min_spikes=args.min_spikes, n_jobs=args.jobs)

    if results.empty:
        print("No valid sessions found.")
    else:
        n_probes = results[['mouse', 'session', 'probe']].drop_duplicates().shape[0]
        print(f"\nDone. {n_probes} probes processed.")
        if args.csv:
            results.to_csv(args.csv, index=False)
            print(f"Raw results saved to {args.csv}")
        plot_summary(results, alpha=args.alpha, save_path=args.save)
