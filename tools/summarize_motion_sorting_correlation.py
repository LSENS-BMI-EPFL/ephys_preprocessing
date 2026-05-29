"""
summarize_motion_sorting_correlation.py

Loads all_recordings_summary.csv produced by run_motion_sorting_correlation.py
and generates dataset-level figures saved to:

  combined_results/motion_sorting_correlation/summary_figures/

Figures
-------
  A  – Paired scatter: fraction non-significant (amplitude ~ motion)
  B  – Paired scatter: fraction non-significant (FR ~ motion)
  C  – r-distribution pooled across all recordings (GOOD and MUA)
  D  – Paired scatter: median presence ratio
  E  – Paired scatter: yield (n units with PR > 0.8)
  F  – Paired scatter: median coverage ratio (GOOD+MUA)
  G  – Fraction of units in over-split quadrant (bar)
  H  – Drift amplitude vs gain in frac_nonsig_amp (GOOD+MUA)
  I  – Drift amplitude vs gain in median presence ratio
  J  – Per-mouse breakdown of key metrics
  K  – Multi-probe consistency (paired imec0 vs imec1)

Usage
-----
  python summarize_motion_sorting_correlation.py [--csv path/to/all_recordings_summary.csv]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu, pearsonr
from itertools import cycle

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")
OUTPUT_BASE    = BASE_DATA_PATH.parent / "combined_results" / "motion_sorting_correlation"
FIG_DIR        = OUTPUT_BASE / "summary_figures"

RUN_COLORS  = {"ks2": "royalblue", "dredge+ks4": "firebrick"}
HIGH_PR     = 0.8
ALPHA_FIG   = 0.01     # for significance annotations on figures

# A colour per mouse (up to ~20)
_MOUSE_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
]


# ---------------------------------------------------------------------------
# Load & pivot
# ---------------------------------------------------------------------------

def load_summary(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Ensure run_label column exists
    if "run_label" not in df.columns:
        df["run_label"] = df["run_name"].apply(
            lambda n: "ks2" if "2" in n.lower().replace("kilosort", "") else "dredge+ks4"
        )
    return df


def pivot_paired(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot so each row is one (mouse, session, probe) with columns for
    both ks2 and dredge+ks4 values side by side.
    Only keeps recordings that have BOTH runs.
    """
    ks2 = df[df["run_label"] == "ks2"].set_index(["mouse", "session", "probe"])
    ks4 = df[df["run_label"] == "dredge+ks4"].set_index(["mouse", "session", "probe"])
    common = ks2.index.intersection(ks4.index)
    ks2 = ks2.loc[common].add_suffix("_ks2")
    ks4 = ks4.loc[common].add_suffix("_ks4")
    paired = pd.concat([ks2, ks4], axis=1).reset_index()
    return paired


def mouse_colors(mice: list) -> dict:
    palette = cycle(_MOUSE_PALETTE)
    dict = {m: next(palette) for m in sorted(set(mice))}
    dict = {m:'k' for m in mice}
    return dict


# ---------------------------------------------------------------------------
# Shared plot helpers
# ---------------------------------------------------------------------------

def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _paired_scatter(
    ax, x: np.ndarray, y: np.ndarray, mouse_labels: list,
    xlabel: str = "ks2", ylabel: str = "dredge+ks4",
    title: str = "",
    lims: tuple | None = None,
    add_wilcoxon: bool = True,
) -> None:
    """Scatter with diagonal, coloured by mouse, optional Wilcoxon annotation."""
    m_colors = mouse_colors(mouse_labels)
    for xi, yi, mi in zip(x, y, mouse_labels):
        ax.scatter(xi, yi, s=18, alpha=0.75, color=m_colors[mi],
                   linewidths=0, zorder=3)

    valid = np.isfinite(x) & np.isfinite(y)
    if lims is None:
        pad  = 0.04 * max(np.nanmax(np.abs(x)), np.nanmax(np.abs(y)), 1)
        vmin = min(np.nanmin(x[valid]), np.nanmin(y[valid])) - pad
        vmax = max(np.nanmax(x[valid]), np.nanmax(y[valid])) + pad
        lims = (vmin, vmax)
    ax.plot(lims, lims, "k--", lw=0.8, alpha=0.4)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(xlabel, fontsize=9); ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    _despine(ax)

    if add_wilcoxon and valid.sum() >= 10:
        diff = y[valid] - x[valid]
        try:
            _, p = wilcoxon(diff)
            direction = "↑" if diff.mean() > 0 else "↓"
            p_str = f"p={p:.2e}" if p < 0.001 else f"p={p:.3f}"
            ax.text(0.04, 0.95, f"Wilcoxon {direction}  {p_str}",
                    transform=ax.transAxes, fontsize=7, va="top")
        except Exception:
            pass


def _add_mouse_legend(fig, mice: list, loc="lower right"):
    m_colors = mouse_colors(mice)
    handles  = [plt.scatter([], [], s=20, color=c, label=m)
                for m, c in m_colors.items()]
    leg = fig.legend(handles=handles, title="Mouse", fontsize=7, title_fontsize=8,
                     loc=loc, frameon=False, ncol=max(1, len(mice) // 8))
    return leg


def _save(fig, name: str) -> None:
    path = FIG_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  Saved {path.name}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure A & B – Paired scatter: fraction non-significant
# ---------------------------------------------------------------------------

def fig_nonsig_paired(paired: pd.DataFrame, metric: str, title_tag: str, fname: str) -> None:
    """
    metric: 'amp' or 'fr'
    Three subplots: GOOD, MUA, GOOD+MUA
    """
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle(f"Fraction non-significant  ({title_tag}  ~ motion)", fontsize=12)

    mice = paired["mouse"].tolist()

    for ax, label_tag, col_suf in zip(
        axes,
        ["GOOD", "MUA", "GOOD+MUA"],
        ["good", "mua", "good_mua"],
    ):
        col = f"frac_nonsig_{metric}_{col_suf}"
        x = paired.get(f"{col}_ks2",  pd.Series(dtype=float)).values.astype(float)
        y = paired.get(f"{col}_ks4",  pd.Series(dtype=float)).values.astype(float)
        _paired_scatter(ax, x, y, mice,
                        xlabel="ks2", ylabel="dredge+ks4",
                        title=label_tag, lims=(-0.02, 1.02))

    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure C – r-distribution pooled
# ---------------------------------------------------------------------------

def fig_r_distributions(df: pd.DataFrame, single_dir: Path, fname: str) -> None:
    """
    Load per-unit CSVs and pool r values across all recordings.
    4 panels: r_amp GOOD | r_amp MUA | r_fr GOOD | r_fr MUA
    """
    all_records = []
    for run_csv in sorted(single_dir.glob("*/*/results_*.csv")):
        parts    = run_csv.parts
        mouse    = parts[-3]
        ses_prob = parts[-2]
        run_name = run_csv.stem.replace("results_", "")
        run_lbl  = "ks2" if "2" in run_name.lower().replace("kilosort", "") else "dredge+ks4"
        try:
            sub = pd.read_csv(run_csv, usecols=["label", "r_amp", "r_fr"])
            sub["run_label"] = run_lbl
            all_records.append(sub)
        except Exception:
            pass

    if not all_records:
        print("  [warn] No per-unit CSVs found — skipping Fig C")
        return

    pool = pd.concat(all_records, ignore_index=True)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    fig.suptitle("Pearson r distribution — pooled across all recordings", fontsize=12)

    panels = [
        ("r_amp", "GOOD", "r(amp ~ motion)  GOOD"),
        ("r_amp", "MUA",  "r(amp ~ motion)  MUA"),
        ("r_fr",  "GOOD", "r(FR ~ motion)  GOOD"),
        ("r_fr",  "MUA",  "r(FR ~ motion)  MUA"),
    ]
    colors = {"ks2": "royalblue", "dredge+ks4": "firebrick"}
    for ax, (col_r, lbl, title) in zip(axes, panels):
        for run_lbl, color in colors.items():
            sub = pool[(pool["run_label"] == run_lbl) & (pool["label"] == lbl)][col_r].dropna()
            if len(sub) == 0:
                continue
            ax.hist(sub, bins=40, range=(-1, 1), density=True,
                    histtype="step", lw=2, color=color, label=run_lbl)
            ax.hist(sub, bins=40, range=(-1, 1), density=True,
                    histtype="stepfilled", alpha=0.12, color=color)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_xlabel("Pearson r", fontsize=9)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(title, fontsize=10)
        leg = ax.legend(fontsize=8); leg.get_frame().set_visible(False)
        _despine(ax)

    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure D – Paired scatter: median presence ratio
# ---------------------------------------------------------------------------

def fig_presence_ratio_paired(paired: pd.DataFrame, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(5, 5))
    fig.suptitle("Median presence ratio", fontsize=12)
    mice = paired["mouse"].tolist()
    x = paired["median_presence_ratio_ks2"].values.astype(float)
    y = paired["median_presence_ratio_ks4"].values.astype(float)
    _paired_scatter(ax, x, y, mice, xlabel="ks2", ylabel="dredge+ks4",
                    title="Median presence ratio", lims=(0, 1.05))
    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure E – Paired scatter: yield (n units with PR > 0.8)
# ---------------------------------------------------------------------------

def fig_yield_paired(paired: pd.DataFrame, fname: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(f"Yield  (n units with presence ratio > {HIGH_PR})", fontsize=12)
    mice = paired["mouse"].tolist()

    x_n = paired["n_high_presence_ks2"].values.astype(float)
    y_n = paired["n_high_presence_ks4"].values.astype(float)
    _paired_scatter(axes[0], x_n, y_n, mice,
                    title="Absolute count", xlabel="ks2", ylabel="dredge+ks4")

    x_f = paired["frac_high_presence_ks2"].values.astype(float)
    y_f = paired["frac_high_presence_ks4"].values.astype(float)
    _paired_scatter(axes[1], x_f, y_f, mice,
                    title="Fraction", xlabel="ks2", ylabel="dredge+ks4",
                    lims=(0, 1.05))

    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure F – Paired scatter: median coverage ratio
# ---------------------------------------------------------------------------

def fig_coverage_paired(paired: pd.DataFrame, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(5, 5))
    fig.suptitle("Median coverage ratio  (GOOD+MUA)", fontsize=12)
    mice = paired["mouse"].tolist()
    x = paired["median_coverage_ratio_ks2"].values.astype(float)
    y = paired["median_coverage_ratio_ks4"].values.astype(float)
    _paired_scatter(ax, x, y, mice, xlabel="ks2", ylabel="dredge+ks4",
                    title="Median coverage ratio", lims=(0, 1.05))
    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure G – Fraction over-split (bar + dots per recording)
# ---------------------------------------------------------------------------

def fig_oversplit(df: pd.DataFrame, paired: pd.DataFrame, fname: str) -> None:
    """Bar of mean frac_oversplit per method, with per-recording dots."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(f"Over-split fraction  (coverage < {HIGH_PR} AND presence < {HIGH_PR},  GOOD+MUA)",
                 fontsize=11)

    colors = {"ks2": "royalblue", "dredge+ks4": "firebrick"}
    mice   = paired["mouse"].tolist()
    m_cols = mouse_colors(mice)

    # Left: bar with SEM + dots
    ax = axes[0]
    for i, (run_lbl, color) in enumerate(colors.items()):
        vals = df[df["run_label"] == run_lbl]["frac_oversplit"].dropna().values
        mean = vals.mean(); sem = vals.std() / np.sqrt(len(vals))
        ax.bar(i, mean, color=color, alpha=0.8, label=run_lbl, width=0.5)
        ax.errorbar(i, mean, yerr=sem, color="k", lw=1.5, capsize=4, zorder=5)
        jitter = np.random.default_rng(i).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(i + jitter, vals, s=8, alpha=0.5, color=color, linewidths=0, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(list(colors.keys()), fontsize=9)
    ax.set_ylabel("Fraction over-split units", fontsize=9)
    ax.set_title("Mean ± SEM  +  per-recording dots", fontsize=9)
    ax.set_ylim(0, 1); _despine(ax)
    leg = ax.legend(fontsize=8); leg.get_frame().set_visible(False)

    # Right: paired scatter (one dot per recording, coloured by mouse)
    ax = axes[1]
    x = paired["frac_oversplit_ks2"].values.astype(float)
    y = paired["frac_oversplit_ks4"].values.astype(float)
    _paired_scatter(ax, x, y, mice,
                    title="Paired (per recording)", xlabel="ks2", ylabel="dredge+ks4",
                    lims=(0, 1.02))

    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure H – Drift amplitude vs gain in frac_nonsig_amp
# ---------------------------------------------------------------------------

def fig_drift_vs_nonsig(paired: pd.DataFrame, fname: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle("Drift amplitude vs gain in fraction non-significant (amplitude)", fontsize=11)

    mice   = paired["mouse"].tolist()
    m_cols = mouse_colors(mice)

    drift = paired["drift_amplitude_um_ks2"].values.astype(float)  # shared per recording

    for ax, col_suf, title in zip(axes, ["good", "mua", "good_mua"], ["GOOD", "MUA", "GOOD+MUA"]):
        gain = (paired[f"frac_nonsig_amp_{col_suf}_ks4"].values.astype(float)
                - paired[f"frac_nonsig_amp_{col_suf}_ks2"].values.astype(float))
        valid = np.isfinite(drift) & np.isfinite(gain)
        for xi, yi, mi in zip(drift[valid], gain[valid], np.array(mice)[valid]):
            ax.scatter(xi, yi, s=18, alpha=0.75, color=m_cols[mi], linewidths=0)
        ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
        if valid.sum() >= 5:
            r, p = pearsonr(drift[valid], gain[valid])
            p_str = f"p={p:.2e}" if p < 0.001 else f"p={p:.3f}"
            ax.text(0.04, 0.95, f"r={r:.2f}  {p_str}",
                    transform=ax.transAxes, fontsize=8, va="top")
            m, b = np.polyfit(drift[valid], gain[valid], 1)
            xfit = np.linspace(drift[valid].min(), drift[valid].max(), 100)
            ax.plot(xfit, m * xfit + b, "k-", lw=1, alpha=0.5)
        ax.set_xlabel("Drift amplitude (µm, ks2 recording)", fontsize=9)
        ax.set_ylabel("Gain in frac_nonsig  (ks4 − ks2)", fontsize=9)
        ax.set_title(title, fontsize=10); _despine(ax)

    _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure I – Drift amplitude vs gain in median presence ratio
# ---------------------------------------------------------------------------

def fig_drift_vs_pr(paired: pd.DataFrame, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5))
    fig.suptitle("Drift amplitude vs gain in median presence ratio", fontsize=12)

    mice   = paired["mouse"].tolist()
    m_cols = mouse_colors(mice)

    drift = paired["drift_amplitude_um_ks2"].values.astype(float)
    gain  = (paired["median_presence_ratio_ks4"].values.astype(float)
             - paired["median_presence_ratio_ks2"].values.astype(float))
    valid = np.isfinite(drift) & np.isfinite(gain)

    for xi, yi, mi in zip(drift[valid], gain[valid], np.array(mice)[valid]):
        ax.scatter(xi, yi, s=18, alpha=0.75, color=m_cols[mi], linewidths=0)
    ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
    if valid.sum() >= 5:
        r, p = pearsonr(drift[valid], gain[valid])
        p_str = f"p={p:.2e}" if p < 0.001 else f"p={p:.3f}"
        ax.text(0.04, 0.95, f"r={r:.2f}  {p_str}",
                transform=ax.transAxes, fontsize=8, va="top")
        m_fit, b = np.polyfit(drift[valid], gain[valid], 1)
        xfit = np.linspace(drift[valid].min(), drift[valid].max(), 100)
        ax.plot(xfit, m_fit * xfit + b, "k-", lw=1, alpha=0.5)
    ax.set_xlabel("Drift amplitude (µm, ks2 recording)", fontsize=9)
    ax.set_ylabel("Gain in median PR  (ks4 − ks2)", fontsize=9)
    _despine(ax); _add_mouse_legend(fig, mice)
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure J – Per-mouse breakdown
# ---------------------------------------------------------------------------

def fig_per_mouse(paired: pd.DataFrame, fname: str) -> None:
    """
    Four metrics, one violin per (mouse, run) pair.
    Metrics: frac_nonsig_amp_good_mua, frac_nonsig_fr_good_mua,
             median_presence_ratio, median_coverage_ratio
    """
    metrics = [
        ("frac_nonsig_amp_good_mua", "Frac non-sig amp\n(GOOD+MUA)"),
        ("frac_nonsig_fr_good_mua",  "Frac non-sig FR\n(GOOD+MUA)"),
        ("median_presence_ratio",    "Median presence ratio"),
        ("median_coverage_ratio",    "Median coverage ratio"),
    ]
    mice   = sorted(paired["mouse"].unique())
    colors = {"ks2": "royalblue", "dredge+ks4": "firebrick"}

    fig, axes = plt.subplots(1, len(metrics), figsize=(16, 5))
    fig.suptitle("Per-mouse breakdown  (ks2 vs dredge+ks4)", fontsize=12)

    for ax, (base_col, ylabel) in zip(axes, metrics):
        positions_ks2 = []
        positions_ks4 = []
        vals_ks2      = []
        vals_ks4      = []

        for j, mouse in enumerate(mice):
            sub = paired[paired["mouse"] == mouse]
            v2  = sub[f"{base_col}_ks2"].dropna().values.astype(float)
            v4  = sub[f"{base_col}_ks4"].dropna().values.astype(float)
            if len(v2):
                positions_ks2.append(j - 0.2); vals_ks2.append(v2)
            if len(v4):
                positions_ks4.append(j + 0.2); vals_ks4.append(v4)

        for positions, vals_list, color, lbl in [
            (positions_ks2, vals_ks2, "royalblue", "ks2"),
            (positions_ks4, vals_ks4, "firebrick", "dredge+ks4"),
        ]:
            if not vals_list:
                continue
            parts = ax.violinplot(vals_list, positions=positions,
                                  showmedians=True, showextrema=False, widths=0.35)
            for body in parts["bodies"]:
                body.set_facecolor(color); body.set_alpha(0.45)
            parts["cmedians"].set_color(color)
            # label hack
            ax.scatter([], [], color=color, s=20, label=lbl)

        ax.set_xticks(range(len(mice))); ax.set_xticklabels(mice, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(ylabel, fontsize=9); ax.set_title(ylabel, fontsize=9)
        leg = ax.legend(fontsize=7); leg.get_frame().set_visible(False)
        _despine(ax)

    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Figure K – Multi-probe consistency
# ---------------------------------------------------------------------------

def fig_multi_probe(paired: pd.DataFrame, fname: str) -> None:
    """
    For sessions with ≥2 probes, compare the per-probe gain
    (ks4 − ks2) across probes of the same session.
    """
    paired["session_key"] = paired["mouse"] + "_" + paired["session"]
    multi = paired.groupby("session_key").filter(lambda g: len(g) >= 2)

    if multi.empty:
        print("  [warn] No sessions with ≥2 probes — skipping Fig K")
        return

    metric_pairs = [
        ("frac_nonsig_amp_good_mua", "Gain frac non-sig amp"),
        ("median_presence_ratio",    "Gain median PR"),
    ]
    sessions = multi["session_key"].unique()
    m_cols   = mouse_colors(multi["mouse"].tolist())

    fig, axes = plt.subplots(1, len(metric_pairs), figsize=(10, 5))
    fig.suptitle("Multi-probe consistency  (gain ks4 − ks2 per probe, paired within session)",
                 fontsize=11)

    for ax, (base_col, ylabel) in zip(axes, metric_pairs):
        for ses in sessions:
            sub  = multi[multi["session_key"] == ses].sort_values("probe")
            gain = (sub[f"{base_col}_ks4"].values.astype(float)
                    - sub[f"{base_col}_ks2"].values.astype(float))
            mouse = sub["mouse"].iloc[0]
            color = m_cols[mouse]
            probes = sub["probe"].tolist()
            # Connect probes of the same session with a line
            for pi, (probe, g) in enumerate(zip(probes, gain)):
                ax.scatter(probe, g, s=20, color=color, alpha=0.8, linewidths=0)
            if len(gain) >= 2:
                ax.plot(probes, gain, color=color, lw=0.8, alpha=0.4)

        ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
        ax.set_xlabel("Probe", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(ylabel, fontsize=10); _despine(ax)

    _add_mouse_legend(fig, multi["mouse"].tolist())
    fig.tight_layout()
    _save(fig, fname)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Dataset-level summary figures")
    parser.add_argument(
        "--csv", type=Path,
        default=OUTPUT_BASE / "all_recordings_summary.csv",
        help="Path to the master summary CSV",
    )
    args = parser.parse_args()

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.csv} …")
    df = load_summary(args.csv)
    print(f"  {len(df)} rows  |  {df['run_label'].value_counts().to_dict()}")

    paired = pivot_paired(df)
    print(f"  {len(paired)} paired recordings (both runs present)")

    single_dir = OUTPUT_BASE / "single"

    print("Generating figures …")

    # A – Paired scatter: fraction non-significant (amplitude)
    fig_nonsig_paired(paired, "amp", "amplitude", "figA_nonsig_amp.png")

    # B – Paired scatter: fraction non-significant (FR)
    fig_nonsig_paired(paired, "fr", "FR", "figB_nonsig_fr.png")

    # C – r-distribution pooled
    fig_r_distributions(df, single_dir, "figC_r_distributions.png")

    # D – Paired scatter: median presence ratio
    fig_presence_ratio_paired(paired, "figD_presence_ratio.png")

    # E – Paired scatter: yield
    fig_yield_paired(paired, "figE_yield.png")

    # F – Paired scatter: median coverage ratio
    fig_coverage_paired(paired, "figF_coverage_ratio.png")

    # G – Fraction over-split
    fig_oversplit(df, paired, "figG_oversplit.png")

    # H – Drift vs gain in non-sig fraction
    fig_drift_vs_nonsig(paired, "figH_drift_vs_nonsig.png")

    # I – Drift vs gain in median PR
    fig_drift_vs_pr(paired, "figI_drift_vs_pr.png")

    # J – Per-mouse breakdown
    fig_per_mouse(paired, "figJ_per_mouse.png")

    # K – Multi-probe consistency
    fig_multi_probe(paired, "figK_multi_probe.png")

    print(f"\nAll figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
