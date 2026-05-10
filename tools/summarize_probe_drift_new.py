"""
Probe drift motion summary across Neuropixels recordings.

Usage
-----
  python summarize_probe_drift.py --mode single             # L1 figures only
  python summarize_probe_drift.py --mode dataset            # L2 figures only
  python summarize_probe_drift.py --mode all                # both (default)
  python summarize_probe_drift.py --mode single --mouse MH028
  python summarize_probe_drift.py --mode single --mouse MH028 --session MH028_20250501_104058

Path conventions
----------------
  Per-recording / per-mouse : L1_OUT/
  Dataset figures / table   : L2_OUT/
"""

from __future__ import annotations

import argparse
import json
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import f as f_dist, pearsonr

from spikeinterface.core.motion import Motion

# ── Config ──────────────────────────────────────────────────────────────────
DATA_ROOT  = Path(r"M:\analysis\Axel_Bisi\data")
COMBINED   = Path(r"M:\analysis\Axel_Bisi\combined_results\motion_estimation")
L1_OUT     = COMBINED / "single"
L2_OUT     = COMBINED
MOUSE_INFO = Path(r"M:\share_internal\Axel_Bisi_Share\dataset_info\joint_probe_insertion_info.xlsx")

THRESHOLDS = (5.0, 10.0)
THR_COLORS = ("darkorange", "firebrick")
THR_LABELS = [f"permissive (>{THRESHOLDS[0]} µm/s)", f"conservative (>{THRESHOLDS[1]} µm/s)"]

CMAP_DISP = "PiYG"
FIG_DPI   = 300
N_WORKERS = 2

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8})

INSERTION_DEPTHS: dict[tuple[str, str], float] = {}

# Array fields stored in each record but excluded from DataFrames / CSV
_ARRAY_FIELDS = ("t", "mean_d", "v", "dep", "mean_abs_disp_by_depth")

# ── Utilities ────────────────────────────────────────────────────────────────

def _despine_top_right(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _despine_top_right_bottom(ax):
    for s in ("top", "bottom", "right"):
        ax.spines[s].set_visible(False)


def _load_insertion_depths() -> dict[tuple[str, str], float]:
    if not MOUSE_INFO.exists():
        warnings.warn(f"Insertion info not found: {MOUSE_INFO}")
        return {}
    df = pd.read_excel(MOUSE_INFO)
    return {(str(r.mouse_name), str(r.probe_id)): float(r.depth)
            for r in df.itertuples(index=False)}


def _load_reward_group() -> dict[str, str]:
    if not MOUSE_INFO.exists():
        return {}
    df = pd.read_excel(MOUSE_INFO)
    return {str(row.mouse_name): str(row.reward_group) for row in df.itertuples(index=False)}


def load_motion_info(folder: Path) -> dict:
    folder = Path(folder)
    info = {}
    with open(folder / "parameters.json") as f:
        info["parameters"] = json.load(f)
    with open(folder / "run_times.json") as f:
        info["run_times"] = json.load(f)
    for name in ("peaks", "peak_locations"):
        p = folder / f"{name}.npy"
        info[name] = np.load(p) if p.exists() else None
    if (folder / "motion").is_dir():
        motion = Motion.load(folder / "motion")
    else:
        required = ["spatial_bins.npy", "temporal_bins.npy", "motion.npy"]
        if all((folder / f).is_file() for f in required):
            warnings.warn("Loading Motion from legacy format.")
            motion = Motion(
                displacement=[np.load(folder / "motion.npy")],
                temporal_bins_s=[np.load(folder / "temporal_bins.npy")],
                spatial_bins_um=np.load(folder / "spatial_bins.npy"),
            )
        else:
            warnings.warn("No motion object found.")
            motion = None
    info["motion"] = motion
    return info


# ── Discovery ────────────────────────────────────────────────────────────────

def discover_recordings(root: Path, mouse_filter=None, session_filter=None):
    for p in sorted(root.glob("*/*/Ephys/catgt_*/*_imec*/dredge_fast/motion")):
        parts   = p.parts
        base    = len(root.parts)
        mouse   = parts[base]
        session = parts[base + 1]
        probe   = int(parts[-3][-1:])
        if mouse_filter   and mouse   != mouse_filter:   continue
        if session_filter and session != session_filter: continue
        yield mouse, session, probe, p


# ── Signal helpers ───────────────────────────────────────────────────────────

def _filter_by_depth(disp: np.ndarray, dep: np.ndarray, mouse: str, probe):
    depth_um = INSERTION_DEPTHS.get((mouse, str(probe)))
    if depth_um is not None:
        mask = dep < depth_um
        return disp[:, mask], dep[mask]
    return disp, dep


def _signals(motion: Motion, mouse: str, probe, seg: int = 0):
    disp = motion.displacement[seg]
    t    = motion.temporal_bins_s[seg]
    dep  = motion.spatial_bins_um
    disp, dep = _filter_by_depth(disp, dep, mouse, probe)
    md  = np.mean(disp, axis=1)
    v   = np.gradient(md, t)
    return t, dep, disp, md, v


def detect_peaks_dual(t, v):
    p_per, _ = find_peaks(np.abs(v), height=THRESHOLDS[0])
    p_con, _ = find_peaks(np.abs(v), height=THRESHOLDS[1])
    return p_per, p_con


def _autocorr(x, max_lag_frac=0.5):
    """Normalized autocorrelation, returns lags [0 .. max_lag_frac*N]."""
    x = x - x.mean()
    n = len(x)
    max_lag = max(1, int(n * max_lag_frac))
    full = np.correlate(x, x, mode="full")[n - 1:][:max_lag]
    return full / full[0] if full[0] != 0 else full


def _fit_exponential(t: np.ndarray, md: np.ndarray):
    """
    Fit md(t) = A * exp(-t / tau) + C.

    Stat-test logic
    ---------------
    We run an F-test between two nested models:
      H0 : md(t) = C                     (constant, 1 parameter — no settling)
      H1 : md(t) = A·exp(-t/τ) + C       (3 parameters — exponential settling)

    F = [(RSS_H0 - RSS_H1) / 2] / [RSS_H1 / (n - 3)]

    Under H0, F ~ F(2, n-3).  A significant p-value (p < 0.05) means the
    exponential component explains variance beyond what a flat baseline does —
    i.e., statistically significant settling is present.

    We additionally require τ > 0 (decay, not growth) and τ < duration
    (the settling completes within the recording window) for τ to be labelled
    'valid'.

    Returns a dict {tau, A, C, r2, p_value, tau_valid}, or None on failure.
    """
    def _exp(t, A, tau, C):
        return A * np.exp(-t / tau) + C

    try:
        popt, _ = curve_fit(
            _exp, t, md,
            p0=[float(md[0] - md[-1]), float(t[-1] / 3), float(md[-1])],
            bounds=([-np.inf, 1e-3, -np.inf], [np.inf, np.inf, np.inf]),
            maxfev=8000,
        )
        A, tau, C = popt
        md_pred = _exp(t, *popt)
        ss_res  = float(np.sum((md - md_pred) ** 2))
        ss_tot  = float(np.sum((md - md.mean()) ** 2))
        r2      = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        n       = len(t)
        F       = ((ss_tot - ss_res) / 2) / (ss_res / (n - 3)) if ss_res > 0 else np.inf
        p_value = float(1 - f_dist.cdf(F, 2, n - 3))
        return dict(tau=float(tau), A=float(A), C=float(C),
                    r2=r2, p_value=p_value,
                    tau_valid=bool((tau > 0) and (tau < t[-1])))
    except Exception:
        return None


def compute_metrics(mouse: str, reward_group: str, session: str, probe, motion: Motion) -> dict:
    t, dep, disp, md, v = _signals(motion, mouse, probe)
    p_per, p_con        = detect_peaks_dual(t, v)
    exp                 = _fit_exponential(t, md)
    return dict(
        mouse=mouse, reward_group=reward_group, session=session, probe=probe,
        duration_s=float(t[-1]),
        max_abs_disp_um=float(np.max(np.abs(md))),
        p2p_um=float(md.max() - md.min()),
        rms_um=float(np.sqrt(np.mean(md ** 2))),
        mean_abs_vel_um_s=float(np.mean(np.abs(v))),
        n_fast_events_conservative=int(len(p_con)),
        n_fast_events_permissive=int(len(p_per)),
        tau=exp["tau"]       if exp else np.nan,
        tau_r2=exp["r2"]     if exp else np.nan,
        tau_pval=exp["p_value"] if exp else np.nan,
        tau_valid=exp["tau_valid"] if exp else False,
        insertion_depth=INSERTION_DEPTHS.get((mouse, str(probe)), np.nan),
        # arrays — kept in memory for L2 plots, excluded from CSV
        t=t, mean_d=md, v=v, dep=dep,
        mean_abs_disp_by_depth=np.mean(np.abs(disp), axis=0),
    )


# ── L1-R: per-recording figure ───────────────────────────────────────────────

def _plot_recording_worker(args):
    mouse, session, probe, dredge_path, out_dir, insertion_depths = args
    global INSERTION_DEPTHS
    INSERTION_DEPTHS = insertion_depths
    try:
        motion = load_motion_info(Path(dredge_path))["motion"]
        if motion is None:
            return f"SKIP {mouse}/{session}/imec{probe}: no motion"
        plot_recording(mouse, session, probe, motion, Path(out_dir))
        return f"OK   {mouse}/{session}/imec{probe}"
    except Exception as e:
        return f"ERR  {mouse}/{session}/imec{probe}: {e}"


def plot_recording(mouse: str, session: str, probe, motion: Motion, out_dir: Path):
    t, dep, disp, md, v = _signals(motion, mouse, probe)
    p_per, p_con = detect_peaks_dual(t, v)

    depth_um = INSERTION_DEPTHS.get((mouse, str(probe)))
    title = f"{mouse}  |  {session}  |  imec{probe}"
    if depth_um is not None:
        title += f"  |  insertion depth {depth_um:.0f} µm"

    fig = plt.figure(figsize=(13, 11))
    gs  = gridspec.GridSpec(4, 1, hspace=0.06, height_ratios=[2.2, 1.3, 1.3, 1.5],
                            top=0.94, bottom=0.10, left=0.08, right=0.93)
    ax_h = fig.add_subplot(gs[0])
    ax_d = fig.add_subplot(gs[1], sharex=ax_h)
    ax_v = fig.add_subplot(gs[2], sharex=ax_h)
    ax_e = fig.add_subplot(gs[3], sharex=ax_h)

    clim = np.max(np.abs(disp)) * 1.05
    im = ax_h.imshow(disp.T, aspect="auto", origin="lower", cmap=CMAP_DISP,
                     extent=(t[0], t[-1], dep[0], dep[-1]),
                     vmin=-clim, vmax=clim, interpolation="nearest")
    divider   = make_axes_locatable(ax_h)
    cax_strip = divider.append_axes("bottom", size="6%", pad=0.05)
    cax_strip.set_axis_off()
    cbax = inset_axes(cax_strip, width="25%", height="80%", loc="center")
    cb   = fig.colorbar(im, cax=cbax, orientation="horizontal")
    cb.set_label("Motion (µm)", fontsize=7); cb.ax.tick_params(labelsize=6)
    if depth_um is not None:
        ax_h.axhline(depth_um, color="k", lw=1.2, ls="--", alpha=0.8)
    ax_h.set_ylabel("Depth (µm)"); ax_h.set_title(title)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    n_d     = disp.shape[1]
    cmap_d  = plt.colormaps["viridis"].resampled(max(n_d, 1))
    leg_idx = set(np.round(np.linspace(0, n_d - 1, min(n_d, 8))).astype(int))
    for i in range(n_d):
        lbl = f"{dep[i]:.0f} µm" if i in leg_idx else None
        ax_d.plot(t, disp[:, i], lw=0.8, alpha=0.45, color=cmap_d(i), label=lbl)
    ax_d.plot(t, md, color="k", lw=1.4, label="mean")
    ax_d.axhline(0, color="gray", lw=0.5, ls="--")
    ax_d.set_ylabel("Motion (µm)")
    ax_d.legend(loc="upper right", fontsize=6, handlelength=1, ncol=2, frameon=False)
    plt.setp(ax_d.get_xticklabels(), visible=False)
    _despine_top_right_bottom(ax_d)

    ax_v.plot(t, v, color="steelblue", lw=0.9)
    ax_v.axhline(0, color="gray", lw=0.5, ls="--")
    ax_v.set_ylabel("Velocity (µm/s)")
    plt.setp(ax_v.get_xticklabels(), visible=False)
    _despine_top_right_bottom(ax_v)

    ax_e.plot(t, v, color="steelblue", lw=0.85, alpha=0.6, label="velocity")
    ax_e2 = ax_e.twinx()
    ax_e2.plot(t, md, color="0.25", lw=1.1, alpha=0.8, label="motion")
    ax_e2.set_ylabel("Motion (µm)", color="0.25")
    _despine_top_right_bottom(ax_e2)

    for peaks_idx, color, label in zip((p_per, p_con), THR_COLORS, THR_LABELS):
        ax_e.scatter(t[peaks_idx], v[peaks_idx], color=color, s=35, zorder=5,
                     label=f"{label}  n={len(peaks_idx)}", linewidths=0)
        for pi in peaks_idx:
            ax_e.axvline(t[pi], color=color, lw=0.7, alpha=0.5)
    for thr, color in zip(THRESHOLDS, THR_COLORS):
        ax_e.axhline( thr, color=color, lw=0.7, ls=":", alpha=0.8)
        ax_e.axhline(-thr, color=color, lw=0.7, ls=":", alpha=0.8)

    ax_e.set_ylabel("Velocity (µm/s)"); ax_e.set_xlabel("Time (s)")
    _despine_top_right_bottom(ax_e)
    h_v, l_v = ax_e.get_legend_handles_labels()
    h_m, l_m = ax_e2.get_legend_handles_labels()
    ax_e.legend(h_v + h_m, l_v + l_m, fontsize=6.5,
                loc="lower center", bbox_to_anchor=(0.5, -0.35),
                ncol=4, handlelength=1, frameon=False)

    fig.align_labels()
    fig.savefig(out_dir / f"{mouse}_{session}_imec{probe}.png",
                dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── L1-M: per-mouse figure ───────────────────────────────────────────────────

def plot_mouse(mouse: str, recs: list[tuple], out_dir: Path):
    n      = len(recs)
    colors = plt.colormaps["tab10"].resampled(max(n, 1))
    fig    = plt.figure(figsize=(13, 10))
    gs     = gridspec.GridSpec(3, 1, hspace=0.35, height_ratios=[2, 2, 2.5],
                               top=0.93, bottom=0.07)
    ax_d = fig.add_subplot(gs[0])
    ax_v = fig.add_subplot(gs[1], sharex=ax_d)
    ax_c = fig.add_subplot(gs[2])

    series, labels = [], []
    for i, (session, probe, motion) in enumerate(recs):
        t, _, _, md, v = _signals(motion, mouse, probe)
        lbl = f"imec{probe} | {session}"
        labels.append(lbl); series.append((t, md))
        c = colors(i)
        ax_d.plot(t, md, color=c, lw=1.2, label=lbl)
        ax_v.plot(t, v,  color=c, lw=0.9, alpha=0.85)

    ax_d.axhline(0, color="gray", lw=0.5, ls="--")
    ax_d.set_ylabel("Mean motion (µm)"); ax_d.set_title(f"{mouse} — all probes")
    ax_d.legend(fontsize=6.5, loc="upper right", ncol=min(3, n), frameon=False)
    plt.setp(ax_d.get_xticklabels(), visible=False)
    _despine_top_right(ax_d)

    ax_v.axhline(0, color="gray", lw=0.5, ls="--")
    ax_v.set_ylabel("Velocity (µm/s)"); ax_v.set_xlabel("Time (s)")
    _despine_top_right(ax_v)

    t_end = min(s[0][-1] for s in series)
    t_com = np.linspace(0, t_end, 500)
    mat   = np.array([np.interp(t_com, s[0], s[1]) for s in series])
    corr  = np.corrcoef(mat)

    im = ax_c.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax_c.set_xticks(range(n)); ax_c.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
    ax_c.set_yticks(range(n)); ax_c.set_yticklabels(labels, fontsize=7)
    for i in range(n):
        for j in range(n):
            ax_c.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=7,
                      color="w" if abs(corr[i,j]) > 0.65 else "k")
    plt.colorbar(im, ax=ax_c, label="Pearson r", pad=0.01, fraction=0.025)
    ax_c.set_title("Inter-probe displacement correlation")
    fig.savefig(out_dir / f"{mouse}_all_probes.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── L2 helpers ────────────────────────────────────────────────────────────────

def _mouse_colors(records):
    mice = sorted(set(r["mouse"] for r in records))
    cmap = plt.colormaps["tab10"].resampled(max(len(mice), 1))
    return {m: cmap(i) for i, m in enumerate(mice)}


def _interp_common(records, n=500):
    t_end = min(r["t"][-1] for r in records)
    t_com = np.linspace(0, t_end, n)
    mat   = np.array([np.interp(t_com, r["t"], r["mean_d"]) for r in records])
    return t_com, mat


def _scalar_df(records):
    return pd.DataFrame([{k: v for k, v in r.items() if k not in _ARRAY_FIELDS}
                         for r in records])


# ── L2 existing figures ───────────────────────────────────────────────────────

def plot_grand_overlay(records, out):
    all_amp = np.concatenate([r["mean_d"] for r in records])
    spacing = np.percentile(np.abs(all_amp), 90) * 0.5

    fig, ax = plt.subplots(figsize=(12, max(5, len(records) // 2 * 0.6)))
    yticks, ylabels = [], []
    for i, r in enumerate(records):
        offset = i * spacing
        ax.plot(r["t"], r["mean_d"] + offset, color="k", lw=1.5, alpha=0.75)
        yticks.append(offset)
        ylabels.append(f"{r['mouse']} | {r['session'][-10:]} | imec{r['probe']}")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean motion (µm)  (+ offset per recording)")
    ax.set_title("Grand overlay — all recordings")
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2A_grand_overlay.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_magnitude_distributions(records, out):
    df = _scalar_df(records)
    metrics = ["max_abs_disp_um", "p2p_um", "rms_um", "mean_abs_vel_um_s"]
    ylabels = ["Max |disp| (µm)", "Peak-to-peak (µm)", "RMS (µm)", "Mean |vel| (µm/s)"]

    fig, axes = plt.subplots(len(metrics), 2, figsize=(10, 12))
    for i, (metric, ylabel) in enumerate(zip(metrics, ylabels)):
        for col, (hue, palette, col_title) in enumerate([
            (None, None, "All mice"),
            ("reward_group", {"R+": "forestgreen", "R-": "crimson"}, "By reward group"),
        ]):
            ax = axes[i, col]
            if hue is None:
                sns.kdeplot(data=df, x=metric, fill=True, color="k",
                            alpha=0.3, lw=1.5, ax=ax)
                sns.rugplot(data=df, x=metric, color="k", height=0.05, ax=ax)
            else:
                sns.kdeplot(data=df, x=metric, hue=hue, hue_order=["R+", "R-"],
                            palette=palette, fill=True, alpha=0.25, lw=1.5,
                            common_norm=False, ax=ax)
                sns.rugplot(data=df, x=metric, hue=hue, hue_order=["R+", "R-"],
                            palette=palette, height=0.05, ax=ax)
            ax.set_ylabel(ylabel); ax.set_xlabel("")
            if i == 0: ax.set_title(col_title)
            _despine_top_right(ax)
    for ax in axes[-1]:
        ax.set_xlabel("Value")
    fig.suptitle("Motion magnitude distributions", y=1.01)
    fig.tight_layout()
    fig.savefig(out / "L2B_magnitude_distributions.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_magnitude_distributions_sessions(records, out):
    df = _scalar_df(records)
    metrics = ["max_abs_disp_um", "p2p_um", "rms_um", "mean_abs_vel_um_s"]
    ylabels = ["Max |disp| (µm)", "Peak-to-peak (µm)", "RMS (µm)", "Mean |vel| (µm/s)"]

    fig, axes = plt.subplots(4, 1, figsize=(8, 14))
    for ax, metric, ylabel in zip(axes, metrics, ylabels):
        sns.violinplot(data=df, x="mouse", y=metric, ax=ax, inner=None,
                       color="0.75", linewidth=0.8)
        sns.stripplot(data=df, x="mouse", y=metric, ax=ax, color="k",
                      size=4, jitter=True, linewidth=0)
        ax.set_ylabel(ylabel); ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=30)
        _despine_top_right(ax)
    axes[-1].set_xlabel("Mouse")
    fig.suptitle("Motion magnitude distributions")
    fig.tight_layout()
    fig.savefig(out / "L2B_magnitude_distributions_sessions.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_drift_correlation(records, out):
    import scipy.stats
    df = _scalar_df(records)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    for ax in axs:
        _despine_top_right(ax)

    ax = axs[0]
    ax.scatter(df["duration_s"], df["max_abs_disp_um"],
               color="k", s=25, alpha=0.8, linewidths=0)
    coeffs = np.polyfit(df["duration_s"], df["max_abs_disp_um"], 1)
    xr = np.array([df["duration_s"].min(), df["duration_s"].max()])
    ax.plot(xr, np.polyval(coeffs, xr), color="k", lw=1.5)
    r, p = scipy.stats.pearsonr(df["duration_s"], df["max_abs_disp_um"])
    ax.annotate(f"r={r:.2f}, p={p:.3f}", xy=(0.05, 0.95), xycoords="axes fraction",
                ha="left", va="top", fontsize=7)
    ax.set_xlabel("Recording duration (s)"); ax.set_ylabel("Max |disp| (µm)")
    ax.set_title("Drift vs. recording duration")

    ax = axs[1]
    palette = {"R+": "forestgreen", "R-": "crimson"}
    sns.scatterplot(data=df, x="duration_s", y="max_abs_disp_um",
                    hue="reward_group", hue_order=["R+", "R-"],
                    palette=palette, ax=ax, s=25, alpha=0.8, edgecolor="none")
    y_pos = 0.95
    for grp, color in palette.items():
        sub = df[df["reward_group"] == grp]
        if len(sub) < 2: continue
        c = np.polyfit(sub["duration_s"], sub["max_abs_disp_um"], 1)
        xr = np.array([sub["duration_s"].min(), sub["duration_s"].max()])
        ax.plot(xr, np.polyval(c, xr), color=color, lw=1.5)
        r, p = scipy.stats.pearsonr(sub["duration_s"], sub["max_abs_disp_um"])
        ax.annotate(f"{grp}: r={r:.2f}, p={p:.3f}", xy=(0.05, y_pos),
                    xycoords="axes fraction", ha="left", va="top",
                    fontsize=7, color=color)
        y_pos -= 0.08
    ax.set_xlabel("Recording duration (s)"); ax.set_ylabel("Max |disp| (µm)")
    ax.set_title("Drift vs. recording duration")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "L2C_drift_vs_duration.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_cross_recording_correlation(records, out):
    _, mat = _interp_common(records, n=1000)
    n    = len(records)
    mice = [r["mouse"] for r in records]

    corr   = np.corrcoef(mat)
    norm_r = plt.Normalize(vmin=-1, vmax=1)
    cmap_r = plt.colormaps["RdBu_r"]
    sm     = plt.cm.ScalarMappable(cmap=cmap_r, norm=norm_r)

    all_x, all_y, all_c     = [], [], []
    same_x, same_y, same_cr = [], [], []
    diff_x, diff_y           = [], []

    for i in range(n):
        for j in range(i + 1, n):
            x, y = mat[i], mat[j]
            col  = cmap_r(norm_r(corr[i, j]))
            all_x.extend(x); all_y.extend(y); all_c.extend([col] * len(x))
            if mice[i] == mice[j]:
                same_x.extend(x); same_y.extend(y); same_cr.extend([col] * len(x))
            else:
                diff_x.extend(x); diff_y.extend(y)

    kw = dict(s=1, alpha=0.1, linewidths=0, rasterized=True)
    for suffix, bxy, bc, oxy, oc, title in [
        ("all",     (all_x, all_y),  all_c,
         None, None, "Pairwise motion — all pairs"),
        ("bymouse", (diff_x, diff_y), ["0.82"] * len(diff_x),
         (same_x, same_y), same_cr, "Pairwise motion — same mouse coloured by r"),
    ]:
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.scatter(*bxy, c=bc, **kw)
        if oxy is not None and oxy[0]:
            ax.scatter(*oxy, c=oc, **kw)
        ax.set_xlabel("Motion recording A (µm)")
        ax.set_ylabel("Motion recording B (µm)")
        ax.set_title(title)
        ax.set_aspect("equal", "datalim")
        _despine_top_right(ax)
        fig.colorbar(sm, ax=ax, label="Pearson r", shrink=0.7, pad=0.02)
        fig.tight_layout()
        fig.savefig(out / f"L2D_pairwise_scatter_{suffix}.png", dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)


def plot_fast_event_summary(records, out):
    df = _scalar_df(records)
    df = df.sort_values(["mouse", "session", "probe"]).reset_index(drop=True)
    lbls = [f"{r.mouse} | {r.session[-15:]} | imec{r.probe}" for r in df.itertuples()]
    y    = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.45)))
    ax.barh(y, df["n_fast_events_permissive"],   color=THR_COLORS[1], alpha=0.55, label=THR_LABELS[1])
    ax.barh(y, df["n_fast_events_conservative"], color=THR_COLORS[0], alpha=0.9,  label=THR_LABELS[0])
    ax.set_yticks(y); ax.set_yticklabels(lbls, fontsize=7)
    ax.set_xlabel("Number of fast events"); ax.set_title("Fast-event count per recording")
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2E_fast_event_summary.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_drift_over_session(records, out):
    by_session = defaultdict(list)
    for r in records:
        by_session[(r["mouse"], r["session"])].append(r)
    for k in by_session:
        by_session[k].sort(key=lambda r: r["probe"])

    fig, ax = plt.subplots(figsize=(7, 5))
    for recs in by_session.values():
        xs = list(range(1, len(recs) + 1))
        ys = [r["mean_abs_vel_um_s"] for r in recs]
        ax.plot(xs, ys, color="0.5", lw=0.9, alpha=0.6)
        ax.scatter(xs, ys, color="k", s=25, zorder=3, linewidths=0)
    ax.set_xlabel("Probe index within session")
    ax.set_ylabel("Mean |velocity| (µm/s)")
    ax.set_title("Probe drift across session")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2C_drift_over_session.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── New L2 figures ────────────────────────────────────────────────────────────

def plot_exponential_decay(records, out):
    """
    Settling time constant τ from fitting md(t) = A·exp(-t/τ) + C per recording.

    Layout
    ------
    Row 0 — τ distribution (all valid fits) | τ by reward group | dataset-wide fit
    Row 1 — 3 example individual fits with annotated τ, R², p
    """
    def _exp(t, A, tau, C): return A * np.exp(-t / tau) + C

    valid = [r for r in records if r["tau_valid"] and not np.isnan(r["tau"])]
    taus  = np.array([r["tau"] for r in valid])
    sig   = np.array([r["tau_pval"] < 0.05 for r in valid])

    fig = plt.figure(figsize=(14, 9))
    gs  = gridspec.GridSpec(2, 3, hspace=0.50, wspace=0.38,
                            top=0.92, bottom=0.08, left=0.07, right=0.97)
    ax_dist  = fig.add_subplot(gs[0, 0])
    ax_group = fig.add_subplot(gs[0, 1])
    ax_grand = fig.add_subplot(gs[0, 2])
    ax_ex    = [fig.add_subplot(gs[1, k]) for k in range(3)]

    # τ KDE — all recordings
    if len(taus):
        sns.kdeplot(taus, fill=True, color="k", alpha=0.3, lw=1.5, ax=ax_dist)
        sns.rugplot(taus, color="k", height=0.05, ax=ax_dist)
        med = np.median(taus)
        ax_dist.axvline(med, color="k", lw=1.2, ls="--")
        ax_dist.annotate(
            f"median τ = {med:.0f} s\n{100*sig.mean():.0f}% significant (p<0.05)",
            xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top", fontsize=7)
    ax_dist.set_xlabel("τ (s)"); ax_dist.set_ylabel("Density")
    ax_dist.set_title("Settling τ — all recordings")
    _despine_top_right(ax_dist)

    # τ by reward group
    df_v = pd.DataFrame([{"tau": r["tau"], "reward_group": r["reward_group"]}
                         for r in valid])
    for grp, color in [("R+", "forestgreen"), ("R-", "crimson")]:
        sub = df_v[df_v["reward_group"] == grp]["tau"].dropna()
        if len(sub) >= 2:
            sns.kdeplot(sub, fill=True, color=color, alpha=0.25, lw=1.5,
                        label=grp, ax=ax_group)
            sns.rugplot(sub, color=color, height=0.05, ax=ax_group)
    ax_group.set_xlabel("τ (s)"); ax_group.set_title("τ by reward group")
    ax_group.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_group)

    # Dataset-wide fit on grand-mean trace
    t_com, mat = _interp_common(records, n=1000)
    grand_mean = mat.mean(axis=0)
    exp_grand  = _fit_exponential(t_com, grand_mean)
    ax_grand.plot(t_com, grand_mean, color="k", lw=1.3, label="grand mean")
    if exp_grand:
        fitted = _exp(t_com, exp_grand["A"], exp_grand["tau"], exp_grand["C"])
        ax_grand.plot(t_com, fitted, color="firebrick", lw=1.5, ls="--",
                      label=f"fit  τ={exp_grand['tau']:.0f} s")
        ax_grand.annotate(
            f"τ = {exp_grand['tau']:.0f} s\nR² = {exp_grand['r2']:.2f}"
            f"\np = {exp_grand['p_value']:.3f}",
            xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
            fontsize=7, color="firebrick")
    ax_grand.set_xlabel("Time (s)"); ax_grand.set_ylabel("Motion (µm)")
    ax_grand.set_title("Dataset-wide exponential fit\n(grand-mean trace)")
    ax_grand.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_grand)

    # 3 example individual fits
    n_ex  = min(3, len(valid))
    picks = [valid[int(i)] for i in np.round(np.linspace(0, len(valid) - 1, n_ex))]
    for k, (ax_k, r) in enumerate(zip(ax_ex, picks)):
        t, md = r["t"], r["mean_d"]
        ax_k.plot(t, md, color="k", lw=1.1, label="data")
        exp = _fit_exponential(t, md)
        if exp:
            fitted = _exp(t, exp["A"], exp["tau"], exp["C"])
            ax_k.plot(t, fitted, color="firebrick", lw=1.4, ls="--",
                      label=f"τ={exp['tau']:.0f} s")
            ax_k.set_title(
                f"{r['mouse']} imec{r['probe']}\n"
                f"τ={exp['tau']:.0f} s  R²={exp['r2']:.2f}  p={exp['p_value']:.3f}",
                fontsize=7)
        else:
            ax_k.set_title(f"{r['mouse']} imec{r['probe']} — fit failed", fontsize=7)
        ax_k.set_xlabel("Time (s)")
        if k == 0: ax_k.set_ylabel("Motion (µm)")
        ax_k.legend(frameon=False, fontsize=6)
        _despine_top_right(ax_k)
    for ax_k in ax_ex[n_ex:]:
        ax_k.set_visible(False)

    fig.suptitle("Exponential settling fit   md(t) = A·exp(−t/τ) + C", fontsize=10)
    fig.savefig(out / "L2G_exponential_decay.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_autocorrelation(records, out):
    """
    Velocity autocorrelation — episodic drift → sharp drop-off;
    slow continuous drift → sustained positive correlation at short lags.

    Layout: 6 example panels (or fewer) + 1 summary (mean ± std).
    """
    dt_list = [float(np.mean(np.diff(r["t"]))) for r in records]
    n_lags  = min(int(0.5 * min(len(r["v"]) for r in records)), 500)

    acorrs, lag_s_list = [], []
    for r, dt in zip(records, dt_list):
        ac = _autocorr(r["v"])[:n_lags]
        acorrs.append(ac)
        lag_s_list.append(np.arange(len(ac)) * dt)

    dt_med  = float(np.median(dt_list))
    lag_com = np.arange(n_lags) * dt_med
    mat_ac  = np.array([np.interp(lag_com, ls, ac)
                        for ls, ac in zip(lag_s_list, acorrs)])
    mu_ac, sd_ac = mat_ac.mean(0), mat_ac.std(0)

    n_ex  = min(6, len(records))
    picks = [int(i) for i in np.round(np.linspace(0, len(records) - 1, n_ex))]

    fig = plt.figure(figsize=(14, 8))
    gs  = gridspec.GridSpec(2, 4, hspace=0.50, wspace=0.38,
                            top=0.93, bottom=0.08, left=0.07, right=0.97)
    ex_axes = [fig.add_subplot(gs[row, col]) for row in range(2) for col in range(3)]
    ax_summ = fig.add_subplot(gs[:, 3])

    for k, idx in enumerate(picks):
        ax = ex_axes[k]
        r  = records[idx]
        ax.plot(lag_s_list[idx], acorrs[idx], color="steelblue", lw=0.9)
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.set_title(f"{r['mouse']} imec{r['probe']}", fontsize=7)
        ax.set_xlabel("Lag (s)"); ax.set_ylabel("Autocorr.")
        _despine_top_right(ax)
    for ax in ex_axes[n_ex:]:
        ax.set_visible(False)

    ax_summ.fill_between(lag_com, mu_ac - sd_ac, mu_ac + sd_ac,
                         color="steelblue", alpha=0.25)
    ax_summ.plot(lag_com, mu_ac, color="steelblue", lw=1.5, label="mean")
    ax_summ.axhline(0, color="gray", lw=0.5, ls="--")
    ax_summ.set_xlabel("Lag (s)"); ax_summ.set_ylabel("Autocorrelation")
    ax_summ.set_title("Summary — all recordings\n(mean ± std)")
    ax_summ.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_summ)

    fig.suptitle("Velocity autocorrelation", fontsize=10)
    fig.savefig(out / "L2H_velocity_autocorrelation.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_depth_gradient(records, out):
    """
    Mean |displacement| as a function of spatial bin depth.
    Individual traces in mouse colours; bold mean ± std across all recordings.
    """
    mc = _mouse_colors(records)

    valid = [r for r in records if len(r["dep"]) >= 2]
    if not valid:
        print("  No depth data — skipping L2I."); return

    dep_min = max(r["dep"][0]  for r in valid)
    dep_max = min(r["dep"][-1] for r in valid)
    if dep_min >= dep_max:
        dep_min = min(r["dep"][0]  for r in valid)
        dep_max = max(r["dep"][-1] for r in valid)
    dep_com = np.linspace(dep_min, dep_max, 200)

    profiles = []
    for r in valid:
        prof = np.interp(dep_com, r["dep"], r["mean_abs_disp_by_depth"])
        profiles.append(prof)

    fig, ax = plt.subplots(figsize=(6, 7))
    for r in valid:
        ax.plot(r["mean_abs_disp_by_depth"], r["dep"],
                color=mc[r["mouse"]], lw=0.7, alpha=0.35)
    if profiles:
        mu = np.array(profiles).mean(0)
        sd = np.array(profiles).std(0)
        ax.fill_betweenx(dep_com, mu - sd, mu + sd, color="k", alpha=0.15)
        ax.plot(mu, dep_com, color="k", lw=2, label="mean ± std")

    ax.set_xlabel("Mean |displacement| (µm)")
    ax.set_ylabel("Depth along probe (µm)")
    ax.set_title("Depth gradient of drift magnitude")
    # Mouse colour legend
    handles = [plt.Line2D([0], [0], color=mc[m], lw=1.5, label=m)
               for m in sorted(mc)]
    handles.append(plt.Line2D([0], [0], color="k", lw=2, label="mean ± std"))
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper right")
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2I_depth_gradient.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_insertion_depth_scatter(records, out):
    """
    Insertion depth vs. drift metrics and settling τ, with Pearson r annotation.
    Tests whether deeper insertions drift more or settle differently.
    """
    df = _scalar_df(records)
    df = df[df["insertion_depth"].notna() & np.isfinite(df["insertion_depth"])].copy()
    if df.empty:
        print("  No insertion depth data — skipping L2J."); return

    pairs = [
        ("p2p_um",            "Peak-to-peak (µm)"),
        ("rms_um",            "RMS (µm)"),
        ("tau",               "Settling τ (s)"),
        ("mean_abs_vel_um_s", "Mean |vel| (µm/s)"),
    ]
    fig, axes = plt.subplots(1, len(pairs), figsize=(14, 4))
    for ax, (metric, ylabel) in zip(axes, pairs):
        sub = df[df[metric].notna() & np.isfinite(df[metric])]
        ax.scatter(sub["insertion_depth"], sub[metric],
                   color="k", s=25, alpha=0.7, linewidths=0)
        if len(sub) >= 3:
            coeffs = np.polyfit(sub["insertion_depth"], sub[metric], 1)
            xr = np.array([sub["insertion_depth"].min(), sub["insertion_depth"].max()])
            ax.plot(xr, np.polyval(coeffs, xr), color="firebrick", lw=1.4, ls="--")
            r, p = pearsonr(sub["insertion_depth"], sub[metric])
            ax.annotate(f"r={r:.2f}, p={p:.3f}", xy=(0.05, 0.95),
                        xycoords="axes fraction", ha="left", va="top",
                        fontsize=7, color="firebrick")
        ax.set_xlabel("Insertion depth (µm)"); ax.set_ylabel(ylabel)
        _despine_top_right(ax)

    fig.suptitle("Drift metrics vs. insertion depth", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "L2J_insertion_depth_scatter.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def save_summary_table(records, out):
    df = _scalar_df(records)
    df = df.sort_values(["mouse", "session", "probe"]).reset_index(drop=True)
    df.to_csv(out / "motion_summary.csv", index=False)

    numeric_cols = df.select_dtypes(include="number").columns
    df[numeric_cols] = df[numeric_cols].round(3)
    fig, ax = plt.subplots(figsize=(18, max(2, len(df) * 0.38 + 0.8)))
    ax.axis("off")
    tbl = ax.table(cellText=df.values.tolist(), colLabels=df.columns.tolist(),
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(6.5)
    tbl.auto_set_column_width(range(len(df.columns)))
    fig.tight_layout()
    fig.savefig(out / "L2F_summary_table.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Probe drift motion summary.")
    parser.add_argument("--mode", choices=["single", "dataset", "all"], default="dataset")
    parser.add_argument("--mouse",   default=None)
    parser.add_argument("--session", default=None)
    parser.add_argument("--workers", type=int, default=N_WORKERS)
    args = parser.parse_args()

    global INSERTION_DEPTHS
    INSERTION_DEPTHS = _load_insertion_depths()
    REWARD_GROUPS    = _load_reward_group()
    L1_OUT.mkdir(parents=True, exist_ok=True)
    L2_OUT.mkdir(parents=True, exist_ok=True)

    all_recs = list(discover_recordings(DATA_ROOT, args.mouse, args.session))
    if not all_recs:
        print("No recordings found."); return

    records:  list[dict]      = []
    by_mouse: dict[str, list] = defaultdict(list)

    print(f"Loading {len(all_recs)} recording(s)...")
    for mouse, session, probe, dredge_path in all_recs:
        reward_group = REWARD_GROUPS.get(mouse, "unknown")
        try:
            motion = load_motion_info(dredge_path)["motion"]
        except FileNotFoundError as e:
            print(f"  SKIP {mouse}/{session}/imec{probe}: {e}"); continue
        if motion is None:
            print(f"  SKIP {mouse}/{session}/imec{probe}: no motion object"); continue
        by_mouse[mouse].append((session, probe, motion))
        records.append(compute_metrics(mouse, reward_group, session, probe, motion))
        print(f"  loaded {mouse} | {session} | imec{probe}")

    if not records:
        print("No valid motion objects found."); return

    if args.mode in ("single", "all"):
        print(f"\nGenerating L1 per-recording figures ({args.workers} workers)...")
        worker_args = [
            (mouse, session, probe, str(dredge_path), str(L1_OUT), INSERTION_DEPTHS)
            for mouse, session, probe, dredge_path in all_recs
        ]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(_plot_recording_worker, a): a for a in worker_args}
            for fut in as_completed(futs):
                print(" ", fut.result())

        print("\nGenerating L1 per-mouse figures...")
        for mouse, recs in by_mouse.items():
            plot_mouse(mouse, recs, L1_OUT)
            print(f"  {mouse}_all_probes.png")

    if args.mode in ("dataset", "all"):
        print("\nGenerating L2 dataset figures...")
        plot_grand_overlay(records, L2_OUT)
        plot_magnitude_distributions(records, L2_OUT)
        plot_magnitude_distributions_sessions(records, L2_OUT)
        plot_drift_correlation(records, L2_OUT)
        plot_cross_recording_correlation(records, L2_OUT)
        plot_fast_event_summary(records, L2_OUT)
        plot_drift_over_session(records, L2_OUT)
        plot_exponential_decay(records, L2_OUT)
        plot_autocorrelation(records, L2_OUT)
        plot_depth_gradient(records, L2_OUT)
        plot_insertion_depth_scatter(records, L2_OUT)
        save_summary_table(records, L2_OUT)

    print(f"\nDone.\n  L1 : {L1_OUT}\n  L2 : {L2_OUT}")


if __name__ == "__main__":
    main()
