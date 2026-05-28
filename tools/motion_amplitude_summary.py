"""
motion_amplitude_summary.py

Run the amplitude–motion correlation analysis across all sessions and probes
found under BASE_DATA_PATH and summarise results in a single figure.

Reuses load_analyzer, compute_correlations, resolve_runs from
motion_amplitude_correlation.py — place both files in the same folder.

Usage
-----
python motion_amplitude_summary.py [--base M:/analysis/Axel_Bisi/data]
                                   [--alpha 0.05] [--min-spikes 50]
                                   [--save summary.png]
"""

import argparse
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import spikeinterface.preprocessing as sp

# Reuse everything from the per-session script
from motion_amplitude_correlation import (
    BASE_DATA_PATH,
    KEEP_LABELS,
    load_analyzer,
    compute_correlations,
)

RUN_COLORS = ['#e69f00', '#56b4e9']


# ---------------------------------------------------------------------------
# Dataset scan
# ---------------------------------------------------------------------------

def iter_sessions(base):
    """
    Yield (mouse, session, probe_name, imec_path, ks_runs) for every valid
    probe found under base. Mirrors resolve_runs logic without raising.
    """
    if not base.exists():
        raise FileNotFoundError(f"Base path not found: {base}")

    for mouse_dir in sorted(base.iterdir()):
        if not mouse_dir.is_dir() or mouse_dir.name.startswith('.'):
            continue
        mouse = mouse_dir.name

        for session_dir in sorted(mouse_dir.iterdir()):
            if not session_dir.is_dir():
                continue
            session = session_dir.name

            for imec in sorted((session_dir / "Ephys").glob("catgt_*/*_imec*")
                               if (session_dir / "Ephys").exists() else []):
                if not imec.is_dir():
                    continue

                motion_path = imec / "dredge_fast" / "motion"
                if not motion_path.is_dir():
                    continue

                imec_part = next(
                    (p for p in reversed(imec.name.split('_')) if p.startswith('imec')),
                    'imec0',
                )

                runs = {}
                for ks_dir in sorted(imec.iterdir()):
                    if not ks_dir.is_dir() or 'kilosort' not in ks_dir.name.lower():
                        continue
                    is_ks2      = ks_dir.name.lower() == "kilosort2"
                    ibl_format  = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
                    metrics     = (ks_dir / "cluster_bc_unitType.tsv" if is_ks2
                                   else ks_dir / "sorter_output" / "cluster_bc_unitType.tsv")
                    if ibl_format.is_dir() and metrics.exists():
                        runs[ks_dir.name] = {
                            'raw_path':    imec,
                            'stream_name': f"{imec_part}.ap",
                            'ibl_format':  ibl_format,
                            'metrics':     metrics,
                            'motion_path': motion_path,
                        }

                if runs:
                    yield mouse, session, imec.name, runs


# ---------------------------------------------------------------------------
# Run across dataset
# ---------------------------------------------------------------------------

def run_dataset(base, alpha, min_spikes):
    """
    Returns a DataFrame with one row per (mouse, session, probe, run_name, label).
    Columns: n_units, n_nonsig, fraction_nonsig.
    Failed sessions are skipped with a warning.
    """
    records = []

    for mouse, session, probe_name, runs in iter_sessions(base):
        tag = f"{mouse} / {session} / {probe_name}"
        print(f"\n{'─'*60}\n{tag}")

        for run_name, paths in runs.items():
            print(f"  {run_name} … ", end="", flush=True)
            try:
                analyzer = load_analyzer(
                    paths['raw_path'], paths['stream_name'],
                    paths['ibl_format'], paths['metrics'],
                )
                motion = sp.load_motion_info(paths['motion_path'])['motion']
                df     = compute_correlations(analyzer, motion, min_spikes=min_spikes)

                for label in ['GOOD', 'MUA']:
                    sub = df[df['label'] == label]
                    n   = len(sub)
                    ns  = int((sub['p_value'] > alpha).sum()) if n else 0
                    records.append({
                        'mouse':          mouse,
                        'session':        session,
                        'probe':          probe_name,
                        'run':            run_name,
                        'label':          label,
                        'n_units':        n,
                        'n_nonsig':       ns,
                        'fraction_nonsig': ns / n if n else np.nan,
                    })
                print(f"{len(df)} units processed")

            except Exception:
                print("FAILED")
                traceback.print_exc()

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Summary figure
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

    # ------------------------------------------------------------------
    # 1. Strip-plot: fraction non-sig per session, one panel per label
    # ------------------------------------------------------------------
    for col, label in enumerate(['GOOD', 'MUA']):
        ax = fig.add_subplot(gs[0, col])
        sub = results[results['label'] == label]

        # One line per session connecting the two run values
        sessions = sub[['mouse', 'session', 'probe']].drop_duplicates()
        for _, row in sessions.iterrows():
            sess_data = sub[
                (sub['mouse'] == row['mouse']) &
                (sub['session'] == row['session']) &
                (sub['probe'] == row['probe'])
            ].set_index('run')
            if len(sess_data) == 2 and run_names[0] in sess_data.index and run_names[1] in sess_data.index:
                y0 = sess_data.loc[run_names[0], 'fraction_nonsig']
                y1 = sess_data.loc[run_names[1], 'fraction_nonsig']
                ax.plot([0, 1], [y0, y1], color='#888', lw=0.8, alpha=0.4, zorder=1)

        # Scatter per run
        for i, run in enumerate(run_names):
            vals = sub[sub['run'] == run]['fraction_nonsig'].dropna()
            jitter = np.random.uniform(-0.04, 0.04, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals,
                       color=colors[run], s=28, alpha=0.75, zorder=2, label=run)
            # Median bar
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

    # ------------------------------------------------------------------
    # 2. Paired difference: KS4 - KS2 fraction non-sig
    # ------------------------------------------------------------------
    ax_diff = fig.add_subplot(gs[0, 2])
    if len(run_names) == 2:
        for label, marker, color in [('GOOD', 'o', '#28aa46'), ('MUA', 's', '#0000dc')]:
            sub = results[results['label'] == label]
            sessions = sub[['mouse', 'session', 'probe']].drop_duplicates()
            diffs = []
            for _, row in sessions.iterrows():
                sess = sub[
                    (sub['mouse'] == row['mouse']) &
                    (sub['session'] == row['session']) &
                    (sub['probe'] == row['probe'])
                ].set_index('run')
                if run_names[0] in sess.index and run_names[1] in sess.index:
                    d = sess.loc[run_names[1], 'fraction_nonsig'] - \
                        sess.loc[run_names[0], 'fraction_nonsig']
                    diffs.append(d)
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
    else:
        ax_diff.text(0.5, 0.5, 'Need exactly 2 runs\nfor paired diff',
                     ha='center', va='center', transform=ax_diff.transAxes,
                     fontsize=10, color='gray')
        ax_diff.set_title('Paired difference')

    # ------------------------------------------------------------------
    # 3. Dataset-level bar: mean fraction per run × label
    # ------------------------------------------------------------------
    ax_bar = fig.add_subplot(gs[1, :2])
    labels   = ['GOOD', 'MUA']
    x        = np.arange(len(labels))
    w        = 0.35
    for i, run in enumerate(run_names):
        means = []
        sems  = []
        for label in labels:
            vals = results[(results['run'] == run) & (results['label'] == label)]['fraction_nonsig'].dropna()
            means.append(vals.mean())
            sems.append(vals.sem())
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

    # ------------------------------------------------------------------
    # 4. Session count table
    # ------------------------------------------------------------------
    ax_tbl = fig.add_subplot(gs[1, 2])
    ax_tbl.axis('off')
    summary_rows = []
    for run in run_names:
        for label in ['GOOD', 'MUA']:
            sub = results[(results['run'] == run) & (results['label'] == label)]
            summary_rows.append([
                run, label,
                f"{len(sub)}",
                f"{sub['n_units'].sum():.0f}",
                f"{sub['fraction_nonsig'].mean():.2f}",
            ])
    tbl = ax_tbl.table(
        cellText=summary_rows,
        colLabels=['Run', 'Label', 'Sessions', 'Units', 'Mean frac'],
        loc='center', cellLoc='center',
    )
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
    parser.add_argument('--save',       default=None)
    parser.add_argument('--csv',        default=None,
                        help='Optional path to save raw per-session results as CSV')
    args = parser.parse_args()

    base = Path(args.base)
    print(f"Scanning {base} …")

    results = run_dataset(base, alpha=args.alpha, min_spikes=args.min_spikes)

    if results.empty:
        print("No valid sessions found.")
    else:
        print(f"\nDone. {len(results)} run×label rows across "
              f"{results[['mouse','session','probe']].drop_duplicates().shape[0]} probes.")
        if args.csv:
            results.to_csv(args.csv, index=False)
            print(f"Raw results saved to {args.csv}")
        plot_summary(results, alpha=args.alpha, save_path=args.save)
