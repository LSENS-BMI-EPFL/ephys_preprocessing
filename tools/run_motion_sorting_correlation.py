"""
run_motion_sorting_correlation.py

Batch runner: iterates over every (mouse, session, probe) under BASE_DATA_PATH,
computes amplitude/FR ~ motion correlations for KS2 and KS4+DREDge, saves:

  combined_results/motion_sorting_correlation/
  ├── single/
  │   └── {mouse}/
  │       └── {session}_{probe}/
  │           ├── figure.png
  │           ├── results_{run_name}.csv   (per-unit DataFrame)
  │           ├── presence_{run_name}.csv
  │           └── done.flag
  └── all_recordings_summary.csv          (one row per mouse/session/probe/run)

Usage
-----
  python run_motion_sorting_correlation.py [--rerun]

  --rerun   reprocess recordings that already have a done.flag
"""

from __future__ import annotations

import argparse
import traceback
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # no display needed in batch mode
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu

import one.alf.io as alfio
import spikeinterface.preprocessing as sp

# ---------------------------------------------------------------------------
# Config — mirrors the GUI script
# ---------------------------------------------------------------------------

BASE_DATA_PATH = Path(r"M:\analysis\Axel_Bisi\data")
OUTPUT_BASE    = BASE_DATA_PATH.parent / "combined_results" / "motion_sorting_correlation"
SINGLE_DIR     = OUTPUT_BASE / "single"

KEEP_LABELS = {"GOOD", "MUA"}
RUN_COLORS  = ["royalblue", "firebrick"]
FR_BIN_SIZE = 2.0   # seconds
ALPHA       = 0.01
MIN_SPIKES  = 50
HIGH_PR     = 0.8   # threshold for yield / bar-chart panels


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _run_label(name: str) -> str:
    return "ks2" if "2" in name.lower().replace("kilosort", "") else "dredge+ks4"


def iter_recordings(base: Path):
    """Yield (mouse, session, probe_str) for every valid imec directory."""
    for mouse_dir in sorted(base.iterdir()):
        if not mouse_dir.is_dir() or mouse_dir.name.startswith("."):
            continue
        for session_dir in sorted(mouse_dir.iterdir()):
            if not session_dir.is_dir():
                continue
            ephys = session_dir / "Ephys"
            if not ephys.is_dir():
                continue
            for imec in sorted(ephys.glob("catgt_*/*_imec*")):
                probe = next(
                    (p for p in reversed(imec.name.split("_")) if p.startswith("imec")),
                    None,
                )
                if probe:
                    yield mouse_dir.name, session_dir.name, probe


def resolve_runs(base: Path, mouse: str, session: str, probe: str) -> dict:
    session_ephys = base / mouse / session / "Ephys"
    imec_dirs = sorted(session_ephys.glob(f"catgt_*/*_{probe}"))
    if not imec_dirs:
        raise FileNotFoundError(f"No imec dir for {probe} in {session_ephys}")
    imec = imec_dirs[0]

    motion_path = imec / "dredge_fast" / "motion"
    if not motion_path.is_dir():
        raise FileNotFoundError(f"Motion folder not found: {motion_path}")

    runs = {}
    for ks_dir in sorted(imec.iterdir()):
        if not ks_dir.is_dir() or "kilosort" not in ks_dir.name.lower():
            continue
        is_ks2      = ks_dir.name.lower() == "kilosort2"
        ibl_format  = imec / "ibl_format" if is_ks2 else ks_dir / "ibl_format"
        metrics     = (ks_dir / "cluster_bc_unitType.tsv" if is_ks2
                       else ks_dir / "sorter_output" / "cluster_bc_unitType.tsv")
        cluster_info = (ks_dir / "cluster_info.tsv" if is_ks2
                        else ks_dir / "sorter_output" / "cluster_info.tsv")
        if ibl_format.is_dir() and metrics.exists():
            runs[ks_dir.name] = {
                "ibl_format":   ibl_format,
                "metrics":      metrics,
                "cluster_info": cluster_info,
                "motion_path":  motion_path,
                "ks_dir":       ks_dir,
                "is_ks2":       is_ks2,
            }
    if not runs:
        raise FileNotFoundError(f"No valid kilosort runs under {imec}")
    return runs


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_spikes(ibl_format: Path, metrics_file: Path, cluster_info_file: Path) -> dict:
    spikes  = alfio.load_object(ibl_format, "spikes")
    metrics = pd.read_csv(metrics_file, sep="\t")
    labels  = metrics["bc_unitType"].values

    ci = pd.read_csv(cluster_info_file, sep="\t")
    cluster_depths = ci["depth"].values

    return {
        "times":          spikes["times"],
        "amps":           spikes["amps"] * 1e6,
        "clusters":       spikes["clusters"].astype(int),
        "labels":         labels,
        "cluster_depths": cluster_depths,
    }


def load_presence_ratio(ks_dir: Path, is_ks2: bool) -> pd.Series | None:
    path = (ks_dir / "qMetrics" / "cluster_presence_ratio.tsv" if is_ks2
            else ks_dir / "sorter_output" / "cluster_presenceRatio.tsv")
    if not path.exists():
        warnings.warn(f"presence ratio not found: {path}")
        return None
    df = pd.read_csv(path, sep="\t")
    ratio_col = next(
        (c for c in df.columns if "ratio" in c.lower() or "presence" in c.lower()),
        df.columns[-1],
    )
    id_col = next(
        (c for c in df.columns if ("id" in c.lower() or "cluster" in c.lower())
         and c != ratio_col), None,
    )
    index = df[id_col].values if id_col else np.arange(len(df))
    return pd.Series(df[ratio_col].values, index=index, name="presence_ratio").dropna()


# ---------------------------------------------------------------------------
# Motion interpolation
# ---------------------------------------------------------------------------

def _interp_motion(motion, times_s: np.ndarray, depth_um: float) -> np.ndarray:
    if not np.isfinite(depth_um):
        return np.full(len(times_s), np.nan)
    t_min, t_max = motion.temporal_bins_s[0][0], motion.temporal_bins_s[0][-1]
    d_min, d_max = motion.spatial_bins_um[0],    motion.spatial_bins_um[-1]
    tc = np.clip(times_s, t_min, t_max)
    dc = float(np.clip(depth_um, d_min, d_max))
    return motion.get_displacement_at_time_and_depth(
        times_s=tc,
        locations_um=np.full(len(tc), dc),
        segment_index=0,
        grid=False,
    ).ravel()


def drift_amplitude_um(motion) -> float:
    """Peak-to-peak of mean displacement across depth, in µm."""
    disp = motion.displacement[0]          # (n_time_bins, n_depth_bins)
    mean_disp = disp.mean(axis=1)
    return float(mean_disp.max() - mean_disp.min())


# ---------------------------------------------------------------------------
# Correlation computation
# ---------------------------------------------------------------------------

def compute_correlations(data: dict, motion, min_spikes: int = MIN_SPIKES) -> pd.DataFrame:
    times          = data["times"]
    amps           = data["amps"]
    clusters       = data["clusters"]
    clusters       = data["clusters"]
    labels         = data["labels"]
    cluster_depths = data["cluster_depths"]

    t_end       = times[-1]
    bin_edges   = np.arange(0.0, t_end + FR_BIN_SIZE, FR_BIN_SIZE)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    records = []
    for uid in np.unique(clusters):
        mask  = clusters == uid
        n     = mask.sum()
        label = labels[uid] if uid < len(labels) else "noise"
        if n < min_spikes:
            continue

        t     = times[mask]
        a     = amps[mask]
        depth = float(cluster_depths[uid]) if uid < len(cluster_depths) else float("nan")
        if not np.isfinite(depth):
            continue

        mot_spikes = _interp_motion(motion, t, depth)
        valid_amp  = np.isfinite(mot_spikes) & np.isfinite(a)
        if valid_amp.sum() >= min_spikes:
            r_amp, p_amp = pearsonr(mot_spikes[valid_amp], a[valid_amp])
        else:
            r_amp, p_amp = np.nan, np.nan

        counts, _ = np.histogram(t, bins=bin_edges)
        fr        = counts / FR_BIN_SIZE
        mot_bins  = _interp_motion(motion, bin_centers, depth)
        valid_fr  = np.isfinite(mot_bins) & np.isfinite(fr)
        if valid_fr.sum() >= 10:
            r_fr, p_fr = pearsonr(mot_bins[valid_fr], fr[valid_fr])
        else:
            r_fr, p_fr = np.nan, np.nan

        coverage_ratio = (t[-1] - t[0]) / t_end if t_end > 0 else np.nan

        records.append({
            "unit_id":        int(uid),
            "label":          label,
            "n_spikes":       int(n),
            "r_amp":          float(r_amp),
            "p_amp":          float(p_amp),
            "r_fr":           float(r_fr),
            "p_fr":           float(p_fr),
            "depth":          depth,
            "coverage_ratio": float(coverage_ratio),
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Scalar summary for one run
# ---------------------------------------------------------------------------

def _scalar_summary(
    mouse: str, session: str, probe: str, run_name: str,
    df: pd.DataFrame, pr: pd.Series | None, motion, alpha: float,
) -> dict:
    """Compute all scalar metrics for one (recording, run) pair."""

    def _frac_ns(sub_df, col_p):
        s = sub_df.dropna(subset=[col_p])
        return float((s[col_p] > alpha).mean()) if len(s) else np.nan

    def _med_abs_r(sub_df, col_r):
        s = sub_df[col_r].dropna()
        return float(s.abs().median()) if len(s) else np.nan

    good    = df[df["label"] == "GOOD"]
    mua     = df[df["label"] == "MUA"]
    good_mua = df[df["label"].isin(KEEP_LABELS)]

    row = {
        "mouse":        mouse,
        "session":      session,
        "probe":        probe,
        "run_name":     run_name,
        "run_label":    _run_label(run_name),
        # drift
        "drift_amplitude_um": drift_amplitude_um(motion),
        "recording_duration_s": float(df["n_spikes"].sum() / max(1, len(df))),  # placeholder
        # unit counts
        "n_units_total":    len(df),
        "n_good":           len(good),
        "n_mua":            len(mua),
        "n_good_mua":       len(good_mua),
        # frac non-sig amplitude
        "frac_nonsig_amp_good":     _frac_ns(good,     "p_amp"),
        "frac_nonsig_amp_mua":      _frac_ns(mua,      "p_amp"),
        "frac_nonsig_amp_good_mua": _frac_ns(good_mua, "p_amp"),
        # frac non-sig FR
        "frac_nonsig_fr_good":      _frac_ns(good,     "p_fr"),
        "frac_nonsig_fr_mua":       _frac_ns(mua,      "p_fr"),
        "frac_nonsig_fr_good_mua":  _frac_ns(good_mua, "p_fr"),
        # median |r| amplitude
        "median_abs_r_amp_good":    _med_abs_r(good,     "r_amp"),
        "median_abs_r_amp_mua":     _med_abs_r(mua,      "r_amp"),
        "median_abs_r_amp_good_mua": _med_abs_r(good_mua, "r_amp"),
        # median |r| FR
        "median_abs_r_fr_good":     _med_abs_r(good,     "r_fr"),
        "median_abs_r_fr_mua":      _med_abs_r(mua,      "r_fr"),
        "median_abs_r_fr_good_mua": _med_abs_r(good_mua, "r_fr"),
        # coverage ratio (GOOD+MUA only)
        "median_coverage_ratio": float(good_mua["coverage_ratio"].median()),
        "n_high_coverage":       int((good_mua["coverage_ratio"] > HIGH_PR).sum()),
        "frac_high_coverage":    float((good_mua["coverage_ratio"] > HIGH_PR).mean())
                                 if len(good_mua) else np.nan,
    }

    # presence ratio (all clusters from TSV, not filtered by label)
    if pr is not None:
        pr_vals = np.asarray(pr.values, dtype=float)
        row.update({
            "n_clusters_pr":        len(pr_vals),
            "median_presence_ratio": float(np.nanmedian(pr_vals)),
            "n_high_presence":       int((pr_vals > HIGH_PR).sum()),
            "frac_high_presence":    float((pr_vals > HIGH_PR).mean()),
        })
        # over-split: join on unit_id to get both metrics
        df_j = df[df["label"].isin(KEEP_LABELS)].copy()
        df_j["presence_ratio"] = df_j["unit_id"].map(pr)
        df_j = df_j.dropna(subset=["coverage_ratio", "presence_ratio"])
        oversplit = (
            (df_j["coverage_ratio"] < HIGH_PR) & (df_j["presence_ratio"] < HIGH_PR)
        )
        row["frac_oversplit"] = float(oversplit.mean()) if len(df_j) else np.nan
    else:
        for k in ("n_clusters_pr", "median_presence_ratio",
                  "n_high_presence", "frac_high_presence", "frac_oversplit"):
            row[k] = np.nan

    return row


# ---------------------------------------------------------------------------
# Plotting (stripped of plt.show; identical logic to GUI script)
# ---------------------------------------------------------------------------

def _sq(ax):
    ax.set_box_aspect(1)
    return ax


def _frac_nonsig(df: pd.DataFrame, col_p: str, alpha: float) -> pd.DataFrame:
    rows = []
    for label in ["GOOD", "MUA", "GOOD+MUA"]:
        sub = (df[df["label"].isin(KEEP_LABELS)] if label == "GOOD+MUA"
               else df[df["label"] == label])
        sub = sub.dropna(subset=[col_p])
        rows.append({
            "label":    label,
            "fraction": float((sub[col_p] > alpha).mean()) if len(sub) else np.nan,
            "n":        len(sub),
        })
    return pd.DataFrame(rows)


def _pick_example_units(df, col_r, col_p, alpha):
    sub = df[df["label"].isin(KEEP_LABELS)].dropna(subset=[col_r, col_p])
    if sub.empty:
        return None, None, None
    abs_r    = sub[col_r].abs()
    uid_corr = int(sub.loc[abs_r.idxmax(), "unit_id"])
    uid_mid  = int(sub.loc[(abs_r - abs_r.median()).abs().idxmin(), "unit_id"])
    nonsig   = sub[sub[col_p] > alpha]
    uid_indep = int(nonsig.loc[nonsig[col_r].abs().idxmin(), "unit_id"]) if not nonsig.empty else None
    return uid_corr, uid_mid, uid_indep


def _scatter_data_for_unit(data, motion, uid):
    mask  = data["clusters"] == uid
    t     = data["times"][mask]
    a     = data["amps"][mask]
    depth = float(data["cluster_depths"][uid])
    mot_spikes = _interp_motion(motion, t, depth)
    t_end       = data["times"][-1]
    bin_edges   = np.arange(0.0, t_end + FR_BIN_SIZE, FR_BIN_SIZE)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts, _   = np.histogram(t, bins=bin_edges)
    fr          = counts / FR_BIN_SIZE
    mot_bins    = _interp_motion(motion, bin_centers, depth)
    return mot_spikes, a, mot_bins, fr


def plot_results(
    results: dict, run_data: dict, run_motions: dict,
    presence: dict, alpha: float,
    session_info: tuple | None = None,
    save_path: Path | None = None,
    show: bool = False,
) -> None:
    run_names  = list(results.keys())
    colors     = RUN_COLORS[: len(run_names)]
    leg_labels = [_run_label(n) for n in run_names]

    ex_run  = run_names[0]
    ex_df   = results[ex_run]
    ex_data = run_data[ex_run]
    ex_mot  = run_motions[ex_run]

    uid_amp_corr, uid_amp_mid, uid_amp_indep = _pick_example_units(ex_df, "r_amp", "p_amp", alpha)
    uid_fr_corr,  uid_fr_mid,  uid_fr_indep  = _pick_example_units(ex_df, "r_fr",  "p_fr",  alpha)

    fig = plt.figure(figsize=(16, 24))
    id_line = " | ".join(session_info) if session_info else ""
    fig.suptitle(f"{id_line} | motion/sorting correlation" , fontsize=11, y=1.0)
    gs = gridspec.GridSpec(6, 4, figure=fig)

    # -- helpers ---
    def _bar_panel(ax, col_p, title):
        fracs  = {n: _frac_nonsig(df, col_p, alpha) for n, df in results.items()}
        labels = ["GOOD", "MUA", "GOOD+MUA"]
        x = np.arange(len(labels)); w = 0.8 / len(run_names)
        for i, (name, color, lbl) in enumerate(zip(run_names, colors, leg_labels)):
            fdf  = fracs[name].set_index("label")
            vals = [fdf.loc[l, "fraction"] if l in fdf.index else np.nan for l in labels]
            ns   = [fdf.loc[l, "n"]        if l in fdf.index else 0       for l in labels]
            bars = ax.bar(x + i * w, vals, width=w * 0.9, color=color, alpha=0.85, label=lbl)
            for bar, n, v in zip(bars, ns, vals):
                if np.isfinite(v):
                    ax.text(bar.get_x() + bar.get_width() / 2, v + 0.02,
                            f"{n}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(x + w * (len(run_names) - 1) / 2); ax.set_xticklabels(labels)
        ax.set_ylim(0, 1.25); ax.set_ylabel("Fraction non-significant  (p > α)")
        ax.set_title(title, fontsize=9)
        leg = ax.legend(fontsize=8); leg.get_frame().set_visible(False)
        ax.axhline(1.0, color="k", lw=0.5, ls="--"); _sq(ax)

    def _hist_panel(ax, col_r, unit_label, title):
        for name, color, lbl in zip(run_names, colors, leg_labels):
            sub = results[name]; sub = sub[sub["label"] == unit_label][col_r].dropna()
            if len(sub) == 0: continue
            ax.hist(sub, bins=30, range=(-1, 1), alpha=0.55,
                    color=color, label=lbl, density=True)
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_xlabel("Pearson r", fontsize=8); ax.set_ylabel("Density", fontsize=8)
        ax.set_title(title, fontsize=9)
        leg = ax.legend(fontsize=7); leg.get_frame().set_visible(False); _sq(ax)

    def _scatter_panel(ax, uid, col_r, col_p, x_vals, y_vals, xlabel, ylabel, ttl):
        if uid is None: ax.set_visible(False); return
        row = ex_df[ex_df["unit_id"] == uid].iloc[0]
        r_val, p_val = row[col_r], row[col_p]
        valid = np.isfinite(x_vals) & np.isfinite(y_vals)
        ax.scatter(x_vals[valid], y_vals[valid], s=4, alpha=0.3, color="k", linewidths=0)
        ax.set_xlabel(xlabel, fontsize=7); ax.set_ylabel(ylabel, fontsize=7)
        ax.tick_params(labelsize=7)
        p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.3f}"
        ax.set_title(f"{ttl}  (u{uid})\nr={r_val:.2f}, p={p_str}", fontsize=7)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    def _draw_distribution_row(row_idx, data_dict, xlabel, threshold, title_prefix):
        active    = [n for n in run_names if data_dict.get(n) is not None and len(data_dict[n]) > 0]
        if not active: return
        ax_hist  = fig.add_subplot(gs[row_idx, 0])
        ax_viol  = fig.add_subplot(gs[row_idx, 1])
        ax_yield = fig.add_subplot(gs[row_idx, 2])
        ax_cnt   = fig.add_subplot(gs[row_idx, 3])
        viol_data = [np.asarray(data_dict[n], dtype=float) for n in active]
        viol_cols = [colors[run_names.index(n)] for n in active]
        viol_lbls = [leg_labels[run_names.index(n)] for n in active]
        positions = np.arange(1, len(active) + 1)
        for name, color, lbl in zip(active, viol_cols, viol_lbls):
            vals = np.asarray(data_dict[name], dtype=float)
            ax_hist.hist(vals, bins=40, range=(0, 1), density=True,
                         histtype="step", lw=1.8, color=color, label=lbl)
            ax_hist.hist(vals, bins=40, range=(0, 1), density=True,
                         histtype="stepfilled", alpha=0.15, color=color)
        ax_hist.set_xlabel(xlabel, fontsize=8); ax_hist.set_ylabel("Density", fontsize=8)
        ax_hist.set_title(f"{title_prefix} distribution", fontsize=9)
        ax_hist.set_xlim(0, 1.05); ax_hist.axvline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
        ax_hist.spines["top"].set_visible(False); ax_hist.spines["right"].set_visible(False)
        leg = ax_hist.legend(fontsize=8, loc="upper left"); leg.get_frame().set_visible(False)
        parts = ax_viol.violinplot(viol_data, positions=positions,
                                   showmedians=False, showextrema=False, widths=0.6)
        for body, color in zip(parts["bodies"], viol_cols):
            body.set_facecolor(color); body.set_alpha(0.35)
            body.set_edgecolor(color); body.set_linewidth(1.2)
        rng = np.random.default_rng(0)
        for pos, vals, color in zip(positions, viol_data, viol_cols):
            q25, med, q75 = np.percentile(vals, [25, 50, 75])
            iqr = q75 - q25
            lo  = max(vals.min(), q25 - 1.5 * iqr); hi = min(vals.max(), q75 + 1.5 * iqr)
            ax_viol.plot([pos, pos], [lo, hi],   color=color, lw=1.2, zorder=3)
            ax_viol.plot([pos, pos], [q25, q75], color=color, lw=5,
                         solid_capstyle="butt", alpha=0.7, zorder=4)
            ax_viol.scatter(pos, med, color="white", s=30, zorder=5,
                            edgecolors=color, linewidths=1.2)
            jitter = rng.uniform(-0.08, 0.08, size=len(vals))
            ax_viol.scatter(pos + jitter, vals, s=3, alpha=0.25,
                            color=color, linewidths=0, zorder=2)
        if len(active) == 2:
            _, p_mw = mannwhitneyu(viol_data[0], viol_data[1], alternative="two-sided")
            p_str   = f"p={p_mw:.2e}" if p_mw < 0.001 else f"p={p_mw:.3f}"
            y_br = 1.07
            ax_viol.plot([1, 1, 2, 2], [y_br - 0.02, y_br, y_br, y_br - 0.02], color="k", lw=0.9)
            ax_viol.text(1.5, y_br + 0.01, f"MW  {p_str}", ha="center", va="bottom", fontsize=7)
        ax_viol.set_xticks(positions); ax_viol.set_xticklabels(viol_lbls, fontsize=8)
        ax_viol.set_ylabel(xlabel, fontsize=8); ax_viol.set_title(f"{title_prefix} by run", fontsize=9)
        ax_viol.set_ylim(-0.05, 1.18); ax_viol.set_xlim(0.4, len(active) + 0.6)
        ax_viol.spines["top"].set_visible(False); ax_viol.spines["right"].set_visible(False)

        thresholds = np.linspace(0, 1, 200)
        for name, color, lbl in zip(active, viol_cols, viol_lbls):
            vals  = np.asarray(data_dict[name], dtype=float)
            fracs = np.array([(vals > thr).sum() / len(vals) for thr in thresholds])
            ax_yield.plot(thresholds, fracs, color=color, lw=2, label=lbl)
        ax_yield.axvline(threshold, color="k", lw=0.8, ls="--", alpha=0.5, label=f"> {threshold}")
        ax_yield.set_xlabel(f"{xlabel} threshold", fontsize=8)
        ax_yield.set_ylabel("Fraction of units", fontsize=8)
        ax_yield.set_title("Yield curve", fontsize=9)
        ax_yield.set_xlim(0, 1); ax_yield.set_ylim(0, 1.05)
        ax_yield.spines["top"].set_visible(False); ax_yield.spines["right"].set_visible(False)
        leg = ax_yield.legend(fontsize=7); leg.get_frame().set_visible(False)
        categories = ["All units", f"> {threshold}"]
        x = np.arange(len(categories)); w = 0.7 / len(active)
        for i, (name, color, lbl) in enumerate(zip(active, viol_cols, viol_lbls)):
            vals   = np.asarray(data_dict[name], dtype=float)
            counts = [len(vals), int((vals > threshold).sum())]
            bars   = ax_cnt.bar(x + i * w, counts, width=w * 0.9, color=color, alpha=0.85, label=lbl)
            for bar, count in zip(bars, counts):
                ax_cnt.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() * 1.02, str(count),
                            ha="center", va="bottom", fontsize=7)
        ax_cnt.set_xticks(x + w * (len(active) - 1) / 2)
        ax_cnt.set_xticklabels(categories, fontsize=8)
        ax_cnt.set_ylabel("# units", fontsize=8)
        ax_cnt.set_title(f"Total vs. high-{title_prefix.lower()} units", fontsize=9)
        ax_cnt.spines["top"].set_visible(False); ax_cnt.spines["right"].set_visible(False)
        leg = ax_cnt.legend(fontsize=7); leg.get_frame().set_visible(False)
        ax_cnt.relim(); ax_cnt.autoscale_view()

    # -- rows --
    _bar_panel(fig.add_subplot(gs[0, :2]), "p_amp",
               "Amplitude ~ motion  |  fraction non-significant")
    _bar_panel(fig.add_subplot(gs[0, 2:]), "p_fr",
               "Firing rate ~ motion  |  fraction non-significant")

    _hist_panel(fig.add_subplot(gs[1, 0]), "r_amp", "GOOD", "r(amp ~ motion)  GOOD")
    _hist_panel(fig.add_subplot(gs[1, 1]), "r_amp", "MUA",  "r(amp ~ motion)  MUA")
    _hist_panel(fig.add_subplot(gs[1, 2]), "r_fr",  "GOOD", "r(FR ~ motion)  GOOD")
    _hist_panel(fig.add_subplot(gs[1, 3]), "r_fr",  "MUA",  "r(FR ~ motion)  MUA")

    gs2 = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[2, :], hspace=1.0, wspace=0.4)
    for col, (uid, ttl) in enumerate(zip(
            [uid_amp_corr, uid_amp_mid, uid_amp_indep], ["correlated", "middle", "independent"])):
        if uid is not None:
            ms, a, _, _ = _scatter_data_for_unit(ex_data, ex_mot, uid)
            _scatter_panel(fig.add_subplot(gs2[0, col]), uid,
                           "r_amp", "p_amp", ms, a, "Motion (µm)", "Amplitude (µV)", f"Amp  {ttl}")
    for col, (uid, ttl) in enumerate(zip(
            [uid_fr_corr, uid_fr_mid, uid_fr_indep], ["correlated", "middle", "independent"])):
        if uid is not None:
            _, _, mb, fr = _scatter_data_for_unit(ex_data, ex_mot, uid)
            _scatter_panel(fig.add_subplot(gs2[1, col]), uid,
                           "r_fr", "p_fr", mb, fr, "Motion (µm)", "Firing rate (Hz)", f"FR  {ttl}")

    valid_runs = [n for n in run_names if presence.get(n) is not None]
    _draw_distribution_row(3, {n: presence[n] for n in valid_runs},
                           xlabel="Presence ratio", threshold=HIGH_PR,
                           title_prefix="Presence ratio")

    cov_data = {
        n: results[n].loc[results[n]["label"].isin(KEEP_LABELS), "coverage_ratio"].dropna().values
        for n in run_names
    }
    _draw_distribution_row(4, cov_data, xlabel="Coverage ratio", threshold=HIGH_PR,
                           title_prefix="Coverage ratio")

    LABEL_COLORS = {"GOOD": "seagreen", "MUA": "darkorange"}
    for col_idx, run_name in enumerate(valid_runs):
        ax_sc = _sq(fig.add_subplot(gs[5, col_idx * 2: col_idx * 2 + 2]))
        pr    = presence[run_name]
        df_r  = results[run_name].copy()
        df_r["presence_ratio"] = df_r["unit_id"].map(pr)
        df_r  = df_r.dropna(subset=["coverage_ratio", "presence_ratio"])
        bg = df_r[~df_r["label"].isin(KEEP_LABELS)]
        if len(bg):
            ax_sc.scatter(bg["coverage_ratio"], bg["presence_ratio"],
                          s=4, alpha=0.15, color="grey", linewidths=0, zorder=1)
        for label, grp in df_r[df_r["label"].isin(KEEP_LABELS)].groupby("label"):
            ax_sc.scatter(grp["coverage_ratio"], grp["presence_ratio"],
                          s=8, alpha=0.5, color=LABEL_COLORS[label],
                          label=label, linewidths=0, zorder=2)
        ax_sc.axvline(HIGH_PR, color="k", lw=0.6, ls="--", alpha=0.4)
        ax_sc.axhline(HIGH_PR, color="k", lw=0.6, ls="--", alpha=0.4)
        for (x, y, txt) in [(0.08, 0.08, "over-split?"),
                            (0.08, 0.90, "dense\nbut short"),
                             (0.84, 0.08, "sparse"),
                            (0.84, 0.90, "stable")]:
            ax_sc.text(x, y, txt, transform=ax_sc.transAxes,
                       fontsize=6, color="#666666", va="bottom" if y > 0.5 else "top")
        ax_sc.set_xlim(-0.02, 1.02); ax_sc.set_ylim(-0.02, 1.08)
        ax_sc.set_xlabel("Coverage ratio  (last−first spike / recording)", fontsize=8)
        ax_sc.set_ylabel("Presence ratio", fontsize=8)
        ax_sc.set_title(f"Over-split diagnostic  –  {leg_labels[run_names.index(run_name)]}", fontsize=9)
        ax_sc.spines["top"].set_visible(False); ax_sc.spines["right"].set_visible(False)
        leg = ax_sc.legend(fontsize=7, markerscale=1.5, loc='upper center'); leg.get_frame().set_visible(False)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Process one recording
# ---------------------------------------------------------------------------

def process_recording(
    mouse: str, session: str, probe: str,
    out_dir: Path, alpha: float, min_spikes: int,
) -> list[dict]:
    """
    Run the full analysis for one (mouse, session, probe).
    Returns a list of scalar summary dicts (one per run).
    """
    runs = resolve_runs(BASE_DATA_PATH, mouse, session, probe)

    results     = {}
    run_data    = {}
    run_motions = {}
    presence    = {}

    for run_name, paths in runs.items():
        print(f"    Loading {run_name} …")
        data   = load_spikes(paths["ibl_format"], paths["metrics"], paths["cluster_info"])
        motion = sp.load_motion_info(paths["motion_path"])["motion"]

        df = compute_correlations(data, motion, min_spikes=min_spikes)
        results[run_name]     = df
        run_data[run_name]    = data
        run_motions[run_name] = motion

        pr = load_presence_ratio(paths["ks_dir"], paths["is_ks2"])
        presence[run_name] = pr

        # Save per-unit CSV
        df.to_csv(out_dir / f"results_{run_name}.csv", index=False)
        if pr is not None:
            pr.reset_index().to_csv(out_dir / f"presence_{run_name}.csv", index=False)

        print(f"      {len(df)} units  |  "
              f"{len(df[df['label'].isin(KEEP_LABELS)])} GOOD/MUA")

    # Save figure
    plot_results(
        results, run_data, run_motions, presence, alpha=alpha,
        session_info=(mouse, session, probe),
        save_path=out_dir / "figure.png",
        show=False,
    )

    # Compute scalar summaries
    summaries = []
    for run_name in runs:
        row = _scalar_summary(
            mouse, session, probe, run_name,
            results[run_name], presence[run_name], run_motions[run_name], alpha,
        )
        summaries.append(row)

    # Save per-recording summary so skipped runs can be recovered later
    pd.DataFrame(summaries).to_csv(out_dir / "summary.csv", index=False)

    return summaries


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Batch motion-sorting correlation")
    parser.add_argument("--rerun", action="store_true",
                        help="Reprocess recordings that already have done.flag")
    parser.add_argument("--alpha", type=float, default=ALPHA,
                        help=f"Significance threshold (default {ALPHA})")
    parser.add_argument("--min-spikes", type=int, default=MIN_SPIKES,
                        help=f"Minimum spikes per unit (default {MIN_SPIKES})")
    args = parser.parse_args()

    SINGLE_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_BASE / "summary_figures").mkdir(parents=True, exist_ok=True)

    all_summaries = []
    recordings    = list(iter_recordings(BASE_DATA_PATH))
    print(f"Found {len(recordings)} recordings to process.")

    #recordings = recordings[::20]
    exclude = ['AB068', 'AB077']
    recordings = [r for r in recordings if r[0] not in exclude]
    #recordings = [r for r in recordings if 'AB126' in r]

    for i, (mouse, session, probe) in enumerate(recordings, 1):
        out_dir = SINGLE_DIR / mouse / f"{session}_{probe}"
        out_dir.mkdir(parents=True, exist_ok=True)
        flag = out_dir / "done.flag"

        if flag.exists() and not args.rerun:
            print(f"[{i}/{len(recordings)}] SKIP {mouse}/{session}/{probe}")
            continue

        print(f"[{i}/{len(recordings)}] {mouse}/{session}/{probe}")
        try:
            summaries = process_recording(
                mouse, session, probe, out_dir,
                alpha=args.alpha, min_spikes=args.min_spikes,
            )
            all_summaries.extend(summaries)
            flag.touch()
        except FileNotFoundError as e:
            print(f"  SKIP — {e}")
        except Exception:
            print(f"  ERROR:")
            traceback.print_exc()

    # Build master CSV: combine freshly computed summaries with any previously
    # saved per-recording summary CSVs from skipped recordings.
    prev_csvs = sorted(SINGLE_DIR.glob("*/*/summary.csv"))
    prev_rows = []
    for csv in prev_csvs:
        try:
            prev_rows.append(pd.read_csv(csv))
        except Exception:
            pass

    all_dfs = ([pd.DataFrame(all_summaries)] if all_summaries else []) + prev_rows
    if all_dfs:
        master  = pd.concat(all_dfs, ignore_index=True).drop_duplicates(
            subset=["mouse", "session", "probe", "run_name"], keep="last"
        )
        out_csv = OUTPUT_BASE / "all_recordings_summary.csv"
        master.to_csv(out_csv, index=False)
        print(f"\nSaved master summary → {out_csv}  ({len(master)} rows)")
    else:
        print("\nNo summaries collected.")


if __name__ == "__main__":
    main()
