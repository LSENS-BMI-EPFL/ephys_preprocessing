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

# TODO: correlate movement with drift...

from __future__ import annotations

import argparse
import json
import warnings
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
from matplotlib.colors import Normalize

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
N_WORKERS = 30

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8})

INSERTION_DEPTHS:  dict[tuple[str, str], float] = {}
TARGET_AREAS:      dict[tuple[str, str], str]   = {}
DAY_OF_RECORDING:  dict[tuple[str, str], int]   = {}

# Array fields excluded from DataFrames / CSV
_ARRAY_FIELDS = ("t", "mean_d", "v", "dep", "mean_abs_disp_by_depth",
                 "whisker_psth", "whisker_t_win", "valve_psth", "valve_t_win")
# ── Utilities ────────────────────────────────────────────────────────────────

TARGET_AREA_CUSTOM_CMAP = {
    'wS1': '#379443', 'wS2': '#51db64', 'A1': '#334b82',
    'PPC': '#5a81db', 'DLS': '#8232ba', 'wM1': '#fab143',
    'wM2': '#ed753e', 'ALM': '#a31010', 'tjM1': '#e3320e',
    'OFC': '#54110c', 'SC':  '#8a6c27',
}

def _despine_top_right(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _despine_top_right_bottom(ax):
    for s in ("top", "bottom", "right"):
        ax.spines[s].set_visible(False)


# ── Metadata loaders ─────────────────────────────────────────────────────────

def _read_metadate_excel(MOUSE_INFO):
    df = pd.read_excel(MOUSE_INFO)
    df = df[df.valid==1]
    return df

def _load_insertion_depths() -> dict[tuple[str, str], float]:
    if not MOUSE_INFO.exists():
        warnings.warn(f"Insertion info not found: {MOUSE_INFO}")
        return {}
    df = _read_metadate_excel(MOUSE_INFO)
    return {(str(r.mouse_name), str(r.probe_id)): float(r.depth)
            for r in df.itertuples(index=False)}


def _load_reward_group() -> dict[str, str]:
    if not MOUSE_INFO.exists():
        return {}
    df = _read_metadate_excel(MOUSE_INFO)
    return {str(row.mouse_name): str(row.reward_group) for row in df.itertuples(index=False)}


def _load_target_area() -> dict[tuple[str, str], str]:
    if not MOUSE_INFO.exists():
        return {}
    df = _read_metadate_excel(MOUSE_INFO)
    if "target_area" not in df.columns:
        warnings.warn("No 'target_area' column in insertion info.")
        return {}
    return {(str(r.mouse_name), str(r.probe_id)): str(r.target_area)
            for r in df.itertuples(index=False)}


def _load_day_of_recording() -> dict[tuple[str, str], int]:
    if not MOUSE_INFO.exists():
        return {}
    df = _read_metadate_excel(MOUSE_INFO)
    # ensure date is a string
    df['date'] = df['date'].astype(str)
    df["date"] = pd.to_datetime(df["date"], format="mixed", dayfirst=True).dt.strftime("%d.%m.%Y")

    if "day_of_recording" not in df.columns:
        warnings.warn("No 'day_of_recording' column in insertion info.")
        return {}
    return {(str(r.mouse_name), str(r.date), str(r.probe_id)): int(r.day_of_recording)
            for r in df.itertuples(index=False)}


# ── Motion loading ───────────────────────────────────────────────────────────

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
    x = x - x.mean()
    n = len(x)
    max_lag = max(1, int(n * max_lag_frac))
    full = np.correlate(x, x, mode="full")[n - 1:][:max_lag]
    return full / full[0] if full[0] != 0 else full


def _fit_exponential(t: np.ndarray, md: np.ndarray, label: str = ""):
    """
    Fit md(t) = A * exp(-t / tau) + C.

    Tries three tau initialisations (fast/medium/slow). Prints failure reason
    when label is provided. Stat-test: F-test vs constant H0.
    """
    def _exp(t, A, tau, C):
        return A * np.exp(-t / tau) + C

    ss_tot = float(np.sum((md - md.mean()) ** 2))
    if ss_tot < 1e-6:
        if label: print(f"    fit skip [{label}]: flat trace")
        return None

    A0, C0, dur = float(md[0] - md[-1]), float(md[-1]), float(t[-1])
    failures = []
    for tau0 in [dur / 10, dur / 3, dur]:
        try:
            popt, _ = curve_fit(
                _exp, t, md,
                p0=[A0, tau0, C0],
                bounds=([-np.inf, 1e-3, -np.inf], [np.inf, np.inf, np.inf]),
                maxfev=12000,
            )
            A, tau, C = popt
            ss_res  = float(np.sum((md - _exp(t, *popt)) ** 2))
            r2      = 1.0 - ss_res / ss_tot
            n       = len(t)
            F       = ((ss_tot - ss_res) / 2) / (ss_res / (n - 3)) if ss_res > 0 else np.inf
            p_value = float(1 - f_dist.cdf(F, 2, n - 3))
            return dict(tau=float(tau), A=float(A), C=float(C),
                        r2=r2, p_value=p_value,
                        tau_valid=bool((tau > 0) and (tau < dur)))
        except Exception as e:
            failures.append(f"tau0={tau0:.0f}s → {e}")
    if label:
        print(f"    fit failed [{label}]: " + "; ".join(failures))
    return None


# ── Whisker stim helpers ──────────────────────────────────────────────────────

def _find_stim_times(dredge_path: Path, filename_stem: str) -> np.ndarray | None:
    """Generic loader for sync event times given a filename stem (no extension)."""
    catgt_dir = dredge_path.parent.parent.parent
    sync_dir  = catgt_dir / "sync_event_times"
    for ext in (".npy", ".txt"):
        p = sync_dir / f"{filename_stem}{ext}"
        if p.exists():
            return np.load(p) if ext == ".npy" else np.loadtxt(p)
    return None


def _find_whisker_times(dredge_path: Path) -> np.ndarray | None:
    return _find_stim_times(dredge_path, "whisker_stim_times")


def _find_valve_times(dredge_path: Path) -> np.ndarray | None:
    return _find_stim_times(dredge_path, "valve_times")


def _compute_psth(t: np.ndarray, md: np.ndarray,
                  stim_times: np.ndarray,
                  pre: float = 2.0, post: float = 5.0,
                  upsample_hz: float = 10.0) -> tuple | None:
    from scipy.interpolate import interp1d

    stim_times = np.asarray(stim_times).flatten()
    stim_times = stim_times[(stim_times > pre) & (stim_times < t[-1] - post)]
    if len(stim_times) == 0:
        return None

    # Cubic interpolant of the full motion trace at upsample_hz
    f_interp = interp1d(t, md, kind="cubic", bounds_error=False, fill_value=np.nan)
    n_pts = int(round((pre + post) * upsample_hz))
    t_win = np.linspace(-pre, post, n_pts)

    trials = []
    for ts in stim_times:
        snippet = f_interp(t_win + float(ts))
        if np.isnan(snippet).any():
            continue
        baseline = snippet[t_win < 0].mean()
        trials.append(snippet - baseline)

    return (np.array(trials), t_win) if trials else None


# ── Metrics ──────────────────────────────────────────────────────────────────

def compute_metrics(mouse: str, reward_group: str, session: str, probe,
                    motion: Motion, dredge_path: Path | None = None) -> dict:
    t, dep, disp, md, v = _signals(motion, mouse, probe)
    p_per, p_con        = detect_peaks_dual(t, v)
    rec_label           = f"{mouse}/{session}/imec{probe}"
    exp                 = _fit_exponential(t, md, label=rec_label)

    whisker_psth = None; whisker_t_win = None
    valve_psth   = None; valve_t_win   = None
    if dredge_path is not None:
        stim_times = _find_whisker_times(Path(dredge_path))
        if stim_times is not None:
            result = _compute_psth(t, md, stim_times)
            if result is not None:
                whisker_psth, whisker_t_win = result
        valve_times = _find_valve_times(Path(dredge_path))
        if valve_times is not None:
            result = _compute_psth(t, md, valve_times)
            if result is not None:
                valve_psth, valve_t_win = result

    session_date = session.split("_")[1]
    YYYY, MM, DD = session_date[:4], session_date[4:6], session_date[6:]
    date = f"{DD}.{MM}.{YYYY}"
    print(date)

    return dict(
        mouse=mouse, reward_group=reward_group, session=session, probe=probe,
        duration_s=float(t[-1]),
        max_abs_disp_um=float(np.max(np.abs(md))),
        p2p_um=float(md.max() - md.min()),
        rms_um=float(np.sqrt(np.mean(md ** 2))),
        mean_abs_vel_um_s=float(np.mean(np.abs(v))),
        n_fast_events_conservative=int(len(p_con)),
        n_fast_events_permissive=int(len(p_per)),
        tau=exp["tau"]          if exp else np.nan,
        tau_r2=exp["r2"]        if exp else np.nan,
        tau_pval=exp["p_value"] if exp else np.nan,
        tau_valid=exp["tau_valid"] if exp else False,
        insertion_depth=INSERTION_DEPTHS.get((mouse, str(probe)), np.nan),
        target_area=TARGET_AREAS.get((mouse, str(probe)), "unknown"),
        day_of_recording=DAY_OF_RECORDING.get((mouse, date, str(probe)), np.nan),
        t=t, mean_d=md, v=v, dep=dep,
        mean_abs_disp_by_depth=np.mean(np.abs(disp), axis=0),
        whisker_psth=whisker_psth,
        whisker_t_win=whisker_t_win,
        valve_psth=valve_psth,
        valve_t_win=valve_t_win,
    )


# ── L1-R: per-recording figure ───────────────────────────────────────────────

def _plot_recording_worker(args):
    mouse, session, probe, dredge_path, out_dir, insertion_depths, target_areas = args
    global INSERTION_DEPTHS, TARGET_AREAS
    INSERTION_DEPTHS = insertion_depths
    TARGET_AREAS     = target_areas
    try:
        info           = load_motion_info(Path(dredge_path))
        motion         = info["motion"]
        if motion is None:
            return f"SKIP {mouse}/{session}/imec{probe}: no motion"
        peaks          = info.get("peaks")
        peak_locations = info.get("peak_locations")
        #plot_recording(mouse, session, probe, motion, Path(out_dir), peaks=peaks, peak_locations=peak_locations)
        return f"OK   {mouse}/{session}/imec{probe}"
    except Exception as e:
        return f"ERR  {mouse}/{session}/imec{probe}: {e}"


def plot_recording(mouse: str, session: str, probe, motion: Motion, out_dir: Path,
                   peaks=None, peak_locations=None):
    t, dep, disp, md, v = _signals(motion, mouse, probe)
    p_per, p_con = detect_peaks_dual(t, v)

    depth_um = INSERTION_DEPTHS.get((mouse, str(probe)))
    title = f"{mouse}  |  {session}  |  imec{probe}"
    if depth_um is not None:
        title += f"  |  insertion depth {depth_um:.0f} µm"

    has_raster = (peaks is not None and peak_locations is not None and len(peaks) > 0)

    if has_raster:
        fig = plt.figure(figsize=(13, 14))
        gs  = gridspec.GridSpec(5, 1, hspace=0.06,
                                height_ratios=[3.0, 1.5, 1.0, 1.0, 1.0],
                                top=0.94, bottom=0.10, left=0.08, right=0.93)
        ax_r = fig.add_subplot(gs[0])
        ax_h = fig.add_subplot(gs[1], sharex=ax_r)
        ax_d = fig.add_subplot(gs[2], sharex=ax_r)
        ax_v = fig.add_subplot(gs[3], sharex=ax_r)
        ax_e = fig.add_subplot(gs[4], sharex=ax_r)
    else:
        fig = plt.figure(figsize=(13, 11))
        gs  = gridspec.GridSpec(4, 1, hspace=0.06,
                                height_ratios=[2.2, 1.3, 1.3, 1.5],
                                top=0.94, bottom=0.10, left=0.08, right=0.93)
        ax_h = fig.add_subplot(gs[0])
        ax_d = fig.add_subplot(gs[1], sharex=ax_h)
        ax_v = fig.add_subplot(gs[2], sharex=ax_h)
        ax_e = fig.add_subplot(gs[3], sharex=ax_h)

    # ── P0 raster ────────────────────────────────────────────────────────────
    if has_raster:
        seg_mask   = peaks["segment_index"] == 0
        pk         = peaks[seg_mask]
        pl         = peak_locations[seg_mask]
        fs_approx  = pk["sample_index"].max() / t[-1] if t[-1] > 0 else 30000
        spike_t    = pk["sample_index"] / fs_approx
        spike_y    = pl["y"]
        amp        = np.abs(pk["amplitude"])
        amp_norm   = (amp - amp.min()) / (np.ptp(amp) + 1e-9)
        if depth_um is not None:
            mask       = spike_y < depth_um
            spike_t    = spike_t[mask]
            spike_y    = spike_y[mask]
            amp_norm   = amp_norm[mask]
        ax_r.scatter(spike_t, spike_y, c=amp_norm, cmap="inferno",
                     s=0.4, alpha=0.3, linewidths=0, rasterized=True)
        ax_r.set_ylabel("Depth (µm)")
        ax_r.set_title(title)
        plt.setp(ax_r.get_xticklabels(), visible=False)
        _despine_top_right_bottom(ax_r)
        title = ""

    # ── P1 heatmap ───────────────────────────────────────────────────────────
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
    ax_h.set_ylabel("Depth (µm)")
    if title:
        ax_h.set_title(title)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # ── P2 traces + mean ─────────────────────────────────────────────────────
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

    # ── P3 velocity ───────────────────────────────────────────────────────────
    ax_v.plot(t, v, color="steelblue", lw=0.9)
    ax_v.axhline(0, color="gray", lw=0.5, ls="--")
    ax_v.set_ylabel("Velocity (µm/s)")
    plt.setp(ax_v.get_xticklabels(), visible=False)
    _despine_top_right_bottom(ax_v)

    # ── P4 fast events ────────────────────────────────────────────────────────
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


# ── L2 figures ───────────────────────────────────────────────────────────────

def plot_grand_overlay(records, out):
    all_amp = np.concatenate([r["mean_d"] for r in records])
    spacing = np.percentile(np.abs(all_amp), 90) * 0.9

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
                            alpha=0.3, lw=1.5, ax=ax, cut=0)
                sns.histplot(data=df, x=metric, bins=50, ax=ax, color="k",
                             alpha=0.3, lw=1.0, stat="density", element="step", fill=False)
                sns.rugplot(data=df, x=metric, color="k", height=0.05, ax=ax)
            else:
                sns.kdeplot(data=df, x=metric, hue=hue, hue_order=["R+", "R-"],
                            palette=palette, fill=True, alpha=0.25, lw=1.5,
                            common_norm=False, ax=ax, cut=0)
                sns.histplot(data=df, x=metric, hue=hue, hue_order=["R+", "R-"],
                             palette=palette, bins=50, ax=ax, alpha=0.3, lw=1.0,
                             stat="density", element="step", fill=False)
                sns.rugplot(data=df, x=metric, hue=hue, hue_order=["R+", "R-"],
                            palette=palette, height=0.05, ax=ax)
            ax.set_ylabel(ylabel); ax.set_xlabel("")
            ax.legend(title="", loc="upper right", frameon=False)
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

    metrics = [
        ("max_abs_disp_um",   "Max |disp| (µm)"),
        ("p2p_um",            "Peak-to-peak (µm)"),
        ("rms_um",            "RMS (µm)"),
        ("mean_abs_vel_um_s", "Mean |vel| (µm/s)"),
        ("tau",               "Settling τ (s)"),
    ]
    palette = {"R+": "forestgreen", "R-": "crimson"}
    df = _scalar_df(records)

    fig, axs = plt.subplots(len(metrics), 2, figsize=(10, 3.5 * len(metrics)))
    for row, (metric, ylabel) in enumerate(metrics):
        needed  = ["duration_s", metric]
        nan_mask = df[needed].isna().any(axis=1)
        for _, r in df[nan_mask].iterrows():
            print(f"  NaN [{metric}]: {r['mouse']} | {r['session']} | imec{r['probe']}")
        sub_all = df.dropna(subset=needed)
        if sub_all.empty:
            for ax in axs[row]: ax.set_visible(False)
            continue

        ax = axs[row, 0]
        ax.scatter(sub_all["duration_s"], sub_all[metric], color="k", s=25, alpha=0.8, linewidths=0)
        if len(sub_all) >= 2:
            coeffs = np.polyfit(sub_all["duration_s"], sub_all[metric], 1)
            xr = np.array([sub_all["duration_s"].min(), sub_all["duration_s"].max()])
            ax.plot(xr, np.polyval(coeffs, xr), color="k", lw=1.5)
            r, p = scipy.stats.pearsonr(sub_all["duration_s"], sub_all[metric])
            ax.annotate(f"r={r:.2f}, p={p:.3f} (n={len(sub_all)})",
                        xy=(0.05, 0.95), xycoords="axes fraction", ha="left", va="top", fontsize=7)
        ax.set_xlabel("Recording duration (s)"); ax.set_ylabel(ylabel)
        if row == 0: ax.set_title("All recordings")
        _despine_top_right(ax)

        ax = axs[row, 1]
        sns.scatterplot(data=sub_all, x="duration_s", y=metric,
                        hue="reward_group", hue_order=["R+", "R-"],
                        palette=palette, ax=ax, s=25, alpha=0.8, edgecolor="none")
        y_pos = 0.95
        for grp, color in palette.items():
            sub = sub_all[sub_all["reward_group"] == grp]
            if len(sub) < 2: continue
            c = np.polyfit(sub["duration_s"], sub[metric], 1)
            xr = np.array([sub["duration_s"].min(), sub["duration_s"].max()])
            ax.plot(xr, np.polyval(c, xr), color=color, lw=1.5)
            r, p = scipy.stats.pearsonr(sub["duration_s"], sub[metric])
            ax.annotate(f"{grp}: r={r:.2f}, p={p:.3f} (n={len(sub)})",
                        xy=(0.05, y_pos), xycoords="axes fraction",
                        ha="left", va="top", fontsize=7, color=color)
            y_pos -= 0.08
        ax.set_xlabel("Recording duration (s)"); ax.set_ylabel(ylabel)
        if row == 0: ax.set_title("By reward group")
        ax.legend(frameon=False, fontsize=7)
        _despine_top_right(ax)

    fig.suptitle("Drift vs. recording duration", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "L2C_drift_vs_duration.png", dpi=FIG_DPI, bbox_inches="tight")
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


def plot_exponential_decay(records, out):
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

    df_v = pd.DataFrame([{"tau": r["tau"], "reward_group": r["reward_group"]} for r in valid])
    for grp, color in [("R+", "forestgreen"), ("R-", "crimson")]:
        sub = df_v[df_v["reward_group"] == grp]["tau"].dropna()
        if len(sub) >= 2:
            sns.kdeplot(sub, fill=True, color=color, alpha=0.25, lw=1.5, label=grp, ax=ax_group)
            sns.rugplot(sub, color=color, height=0.05, ax=ax_group)
    ax_group.set_xlabel("τ (s)"); ax_group.set_title("τ by reward group")
    ax_group.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_group)

    t_com, mat = _interp_common(records, n=1000)
    grand_mean = mat.mean(axis=0)
    exp_grand  = _fit_exponential(t_com, grand_mean)
    ax_grand.plot(t_com, grand_mean, color="k", lw=1.3, label="grand mean")
    if exp_grand:
        fitted = _exp(t_com, exp_grand["A"], exp_grand["tau"], exp_grand["C"])
        ax_grand.plot(t_com, fitted, color="firebrick", lw=1.5, ls="--",
                      label=f"fit  τ={exp_grand['tau']:.0f} s")
        ax_grand.annotate(
            f"τ = {exp_grand['tau']:.0f} s\nR² = {exp_grand['r2']:.2f}\np = {exp_grand['p_value']:.3f}",
            xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
            fontsize=7, color="firebrick")
    ax_grand.set_xlabel("Time (s)"); ax_grand.set_ylabel("Motion (µm)")
    ax_grand.set_title("Dataset-wide exponential fit\n(grand-mean trace)")
    ax_grand.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_grand)

    n_ex  = min(3, len(valid))
    picks = [valid[int(i)] for i in np.round(np.linspace(0, len(valid) - 1, n_ex))]
    for k, (ax_k, r) in enumerate(zip(ax_ex, picks)):
        t, md = r["t"], r["mean_d"]
        ax_k.plot(t, md, color="k", lw=1.1, label="data")
        exp = _fit_exponential(t, md)
        if exp:
            fitted = _exp(t, exp["A"], exp["tau"], exp["C"])
            ax_k.plot(t, fitted, color="firebrick", lw=1.4, ls="--", label=f"τ={exp['tau']:.0f} s")
            ax_k.set_title(
                f"{r['mouse']} imec{r['probe']}\nτ={exp['tau']:.0f} s  R²={exp['r2']:.2f}  p={exp['p_value']:.3f}",
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
    dt_list = [float(np.mean(np.diff(r["t"]))) for r in records]
    n_lags  = min(int(0.5 * min(len(r["mean_d"]) for r in records)), 500)

    acorrs, lag_s_list = [], []
    for r, dt in zip(records, dt_list):
        ac = _autocorr(r["mean_d"])[:n_lags]
        acorrs.append(ac)
        lag_s_list.append(np.arange(len(ac)) * dt)

    dt_med  = float(np.median(dt_list))
    lag_com = np.arange(n_lags) * dt_med
    mat_ac  = np.array([np.interp(lag_com, ls, ac) for ls, ac in zip(lag_s_list, acorrs)])
    mu_ac, sd_ac = mat_ac.mean(0), mat_ac.std(0)

    n_ex  = min(6, len(records))
    picks = [int(i) for i in np.round(np.linspace(0, len(records) - 1, n_ex))]

    fig = plt.figure(figsize=(14, 8))
    gs  = gridspec.GridSpec(2, 4, hspace=0.50, wspace=0.38,
                            top=0.93, bottom=0.08, left=0.07, right=0.97)
    ex_axes = [fig.add_subplot(gs[row, col]) for row in range(2) for col in range(3)]
    ax_summ = fig.add_subplot(gs[:, 3])

    for k, idx in enumerate(picks):
        ax = ex_axes[k]; r = records[idx]
        ax.plot(lag_s_list[idx], acorrs[idx], color="steelblue", lw=0.9)
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.set_title(f"{r['mouse']} imec{r['probe']}", fontsize=7)
        ax.set_xlabel("Lag (s)"); ax.set_ylabel("Autocorr.")
        _despine_top_right(ax)
    for ax in ex_axes[n_ex:]:
        ax.set_visible(False)

    ax_summ.fill_between(lag_com, mu_ac - sd_ac, mu_ac + sd_ac, color="steelblue", alpha=0.25)
    ax_summ.plot(lag_com, mu_ac, color="steelblue", lw=1.5, label="mean")
    ax_summ.axhline(0, color="gray", lw=0.5, ls="--")
    ax_summ.set_xlabel("Lag (s)"); ax_summ.set_ylabel("Autocorrelation")
    ax_summ.set_title("Summary — all recordings\n(mean ± std)")
    ax_summ.legend(frameon=False, fontsize=7)
    _despine_top_right(ax_summ)

    fig.suptitle("Motion autocorrelation", fontsize=10)
    fig.savefig(out / "L2H_motion_autocorrelation.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_depth_gradient(records, out):
    valid = [r for r in records if len(r["dep"]) >= 2]
    if not valid:
        print("  No depth data — skipping L2I."); return

    all_dep = np.concatenate([r["dep"] for r in valid])
    all_val = np.concatenate([r["mean_abs_disp_by_depth"] for r in valid])

    dep_min, dep_max = all_dep.min(), all_dep.max()
    edges   = np.linspace(dep_min, dep_max, 101)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mu, sd  = np.full(100, np.nan), np.full(100, np.nan)
    for k in range(100):
        mask = (all_dep >= edges[k]) & (all_dep < edges[k + 1])
        if mask.sum() >= 2:
            mu[k] = all_val[mask].mean()
            sd[k] = all_val[mask].std()

    fig, ax = plt.subplots(figsize=(6, 7))
    for r in valid:
        ax.plot(r["mean_abs_disp_by_depth"], r["dep"], color="k", lw=0.7, alpha=0.35)

    valid_bins = ~np.isnan(mu)
    ax.fill_betweenx(centers[valid_bins], (mu - sd)[valid_bins], (mu + sd)[valid_bins],
                     color="k", alpha=0.15)
    ax.plot(mu[valid_bins], centers[valid_bins], color="k", lw=2, label="mean ± std")
    ax.set_xlabel("Mean |displacement| (µm)")
    ax.set_ylabel("Depth along probe (µm)")
    ax.set_title("Depth gradient of drift magnitude")
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2I_depth_gradient.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_depth_gradient_absolute(records, out):
    valid = [r for r in records if len(r["dep"]) >= 2 and not np.isnan(r["insertion_depth"])]
    if not valid:
        print("  No absolute depth data — skipping L2I_abs."); return

    all_dep = np.concatenate([(r["insertion_depth"] - r["dep"]) for r in valid])
    all_val = np.concatenate([r["mean_abs_disp_by_depth"] for r in valid])

    dep_min, dep_max = all_dep.min(), all_dep.max()
    edges   = np.linspace(dep_min, dep_max, 51)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mu = np.full(50, np.nan); sd = np.full(50, np.nan)
    for k in range(50):
        mask = (all_dep >= edges[k]) & (all_dep < edges[k + 1])
        if mask.sum() >= 2:
            mu[k] = all_val[mask].mean(); sd[k] = all_val[mask].std()
    vb = ~np.isnan(mu)

    g = sns.JointGrid(x=all_val, y=all_dep, height=8, ratio=4)
    sns.kdeplot(x=all_val, y=all_dep, ax=g.ax_joint, fill=True, cmap="Greys", thresh=0.02, levels=12)
    g.ax_joint.fill_betweenx(centers[vb], (mu - sd)[vb], (mu + sd)[vb], color="firebrick", alpha=0.2)
    g.ax_joint.plot(mu[vb], centers[vb], color="firebrick", lw=2, label="mean ± std")
    g.ax_joint.legend(frameon=False, fontsize=7, loc="lower right")
    g.ax_joint.invert_yaxis()
    g.ax_joint.set_xlabel("Mean |displacement| (µm)")
    g.ax_joint.set_ylabel("Depth below brain surface (µm)")
    g.ax_joint.set_title("Drift magnitude vs. absolute depth", pad=12)
    _despine_top_right(g.ax_joint)

    sns.histplot(x=all_val, ax=g.ax_marg_x, color="k", alpha=0.4, bins=40,
                 stat="density", element="step", fill=True)
    sns.kdeplot(x=all_val, ax=g.ax_marg_x, color="k", lw=1.5)
    sns.rugplot(x=all_val, ax=g.ax_marg_x, color="k", alpha=0.3, height=0.1)
    _despine_top_right(g.ax_marg_x)

    sns.histplot(y=all_dep, ax=g.ax_marg_y, color="k", alpha=0.4, bins=40,
                 stat="density", element="step", fill=True)
    sns.kdeplot(y=all_dep, ax=g.ax_marg_y, color="k", lw=1.5)
    sns.rugplot(y=all_dep, ax=g.ax_marg_y, color="k", alpha=0.3, height=0.1)
    g.ax_marg_y.invert_yaxis()
    _despine_top_right(g.ax_marg_y)

    g.figure.savefig(out / "L2I_depth_gradient_absolute.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(g.figure)


def plot_depth_gradient_barh(records, out):
    valid = [r for r in records if len(r["dep"]) >= 2 and not np.isnan(r["insertion_depth"])]
    if not valid:
        print("  No depth data — skipping depth gradient bar chart."); return

    all_dep = np.concatenate([(r["insertion_depth"] - r["dep"]) for r in valid])
    all_val = np.concatenate([r["mean_abs_disp_by_depth"] for r in valid])

    dep_min, dep_max = all_dep.min(), all_dep.max()
    n_bins  = 40
    edges   = np.linspace(dep_min, dep_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mu = np.full(n_bins, np.nan); sd = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)
    for k in range(n_bins):
        mask = (all_dep >= edges[k]) & (all_dep < edges[k + 1])
        counts[k] = mask.sum()
        if counts[k] >= 2:
            mu[k] = all_val[mask].mean(); sd[k] = all_val[mask].std()

    vb   = ~np.isnan(mu)
    y    = centers[vb]; mu_v = mu[vb]; sd_v = sd[vb]

    fig, ax = plt.subplots(figsize=(6, max(5, n_bins * 0.2)))
    ax.barh(y, mu_v, height=(dep_max - dep_min) / n_bins * 0.8,
            color="k", alpha=0.85, linewidth=0)
    ax.errorbar(mu_v, y, xerr=sd_v, fmt="none", color="k", lw=0.8, capsize=2, alpha=0.6)
    grand_mean = np.nanmean(mu_v)
    ax.axvline(grand_mean, color="firebrick", lw=1.2, ls="--",
               label=f"grand mean = {grand_mean:.2f} µm")
    ax.invert_yaxis()
    ax.set_xlabel("Mean |displacement| (µm)")
    ax.set_ylabel("Depth below brain surface (µm)")
    ax.set_title("Drift magnitude by depth bin")
    ax.legend(frameon=False, fontsize=7)
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2I_depth_gradient_barh.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_insertion_depth_scatter(records, out):
    df = _scalar_df(records)
    df = df[df["insertion_depth"].notna() & np.isfinite(df["insertion_depth"])].copy()
    if df.empty:
        print("  No insertion depth data — skipping L2J."); return

    pairs = [
        ("max_abs_disp_um", "Max |disp| (µm)"),
        ("p2p_um",          "Peak-to-peak (µm)"),
        ("rms_um",          "RMS (µm)"),
        ("tau",             "Settling τ (s)"),
        ("mean_abs_vel_um_s", "Mean |vel| (µm/s)"),
    ]
    fig, axes = plt.subplots(1, len(pairs), figsize=(14, 4))
    for ax, (metric, ylabel) in zip(axes, pairs):
        sub = df[df[metric].notna() & np.isfinite(df[metric])]
        ax.scatter(sub["insertion_depth"], sub[metric], color="k", s=25, alpha=0.7, linewidths=0)
        if len(sub) >= 3:
            coeffs = np.polyfit(sub["insertion_depth"], sub[metric], 1)
            xr = np.array([sub["insertion_depth"].min(), sub["insertion_depth"].max()])
            ax.plot(xr, np.polyval(coeffs, xr), color="firebrick", lw=1.4, ls="--")
            r, p = pearsonr(sub["insertion_depth"], sub[metric])
            ax.annotate(f"r={r:.2f}, p={p:.3f} (n={len(sub)})", xy=(0.05, 0.95),
                        xycoords="axes fraction", ha="left", va="top",
                        fontsize=7, color="firebrick")
        ax.set_xlabel("Insertion depth (µm)"); ax.set_ylabel(ylabel)
        _despine_top_right(ax)

    fig.suptitle("Drift metrics vs. insertion depth", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "L2J_insertion_depth_scatter.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_motion_by_target(records, out):
    df = _scalar_df(records)
    df = df[df["target_area"].notna() & (df["target_area"] != "unknown")].copy()
    if df.empty:
        print("  No target area data — skipping L2K."); return

    metrics = [
        ("max_abs_disp_um",          "Max |disp| (µm)"),
        ("p2p_um",                   "Peak-to-peak (µm)"),
        ("rms_um",                   "RMS (µm)"),
        ("mean_abs_vel_um_s",        "Mean |vel| (µm/s)"),
        ("n_fast_events_permissive", "Fast events (permissive)"),
        ("tau",                      "Settling τ (s)"),
    ]
    metrics = [(m, lbl) for m, lbl in metrics if m in df.columns]

    present = df["target_area"].unique().tolist()
    ordered = [t for t in TARGET_AREA_CUSTOM_CMAP if t in present]
    ordered += [t for t in present if t not in ordered]
    palette = {t: TARGET_AREA_CUSTOM_CMAP.get(t, "#888888") for t in ordered}

    n_metrics = len(metrics)
    fig, axes = plt.subplots(n_metrics, 1,
                             figsize=(max(8, len(ordered) * 0.7), n_metrics * 2.8))
    if n_metrics == 1:
        axes = [axes]

    for ax, (metric, ylabel) in zip(axes, metrics):
        sub = df[df[metric].notna() & np.isfinite(df[metric])]
        sns.stripplot(data=sub, x="target_area", y=metric, order=ordered,
                      palette=palette, ax=ax, size=4, alpha=0.5,
                      jitter=True, linewidth=0, zorder=2)
        for k, target in enumerate(ordered):
            vals = sub[sub["target_area"] == target][metric].dropna()
            if len(vals) < 1: continue
            mu_v = vals.mean(); sd_v = vals.std() if len(vals) > 1 else 0.0
            ax.plot([k - 0.1, k + 0.1], [mu_v, mu_v], color="k", lw=2, zorder=3)
            ax.errorbar(k, mu_v, yerr=sd_v, fmt="none", color="k", lw=1, capsize=0, zorder=3)
        ax.set_ylabel(ylabel); ax.set_xlabel("")
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        _despine_top_right(ax)

    axes[-1].set_xlabel("Target area")
    axes[-1].tick_params(axis="x", rotation=30)
    fig.suptitle("Motion metrics by insertion target", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "L2K_motion_by_target.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_learning_vs_expert(records, out):
    """Compare drift metrics between learning (day 0) and expert (day > 0)."""
    df = _scalar_df(records)
    if "day_of_recording" not in df.columns or df["day_of_recording"].isna().all():
        print("  No day_of_recording data — skipping L2L."); return

    df["stage"] = df["day_of_recording"].apply(
        lambda d: "learning" if pd.notna(d) and int(d) == 0 else "expert"
    )
    # Count number of insertions for each stage:
    print(df.stage.unique(), df["stage"].nunique(), df.stage.value_counts())

    metrics = [
        ("max_abs_disp_um",   "Max |disp| (µm)"),
        ("p2p_um",            "Peak-to-peak (µm)"),
        ("rms_um",            "RMS (µm)"),
        ("mean_abs_vel_um_s", "Mean |vel| (µm/s)"),
        ("n_fast_events_permissive", "Fast events"),
        ("tau",               "Settling τ (s)"),
    ]
    metrics = [(m, lbl) for m, lbl in metrics if m in df.columns]
    palette  = {"learning": "steelblue", "expert": "firebrick"}

    from scipy.stats import mannwhitneyu
    fig, axes = plt.subplots(1, len(metrics), figsize=(3 * len(metrics), 4))
    if len(metrics) == 1: axes = [axes]

    for ax, (metric, ylabel) in zip(axes, metrics):
        sub = df.dropna(subset=[metric])
        sns.boxplot(data=sub, x="stage", y=metric, order=["learning", "expert"],
                    palette=palette, ax=ax, inner=None, linewidth=0.8, alpha=0.5)
        sns.stripplot(data=sub, x="stage", y=metric, order=["learning", "expert"],
                      palette=palette, ax=ax, size=3, jitter=True, linewidth=0, alpha=0.84)
        a = sub[sub["stage"] == "learning"][metric].dropna()
        b = sub[sub["stage"] == "expert"][metric].dropna()
        if len(a) >= 3 and len(b) >= 3:
            _, p = mannwhitneyu(a, b, alternative="two-sided")
            ax.set_title(f"p={p:.3f}", fontsize=7)
        ax.set_ylabel(ylabel); ax.set_xlabel("")
        _despine_top_right(ax)

    fig.suptitle("Learning (day 0) vs. expert sessions")
    fig.tight_layout()
    fig.savefig(out / "L2L_learning_vs_expert.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_early_drift(records, out, window_s=1800):
    """Drift time course over first window_s seconds. Panel A: overlay. Panel B: mean ± std."""
    clipped = []
    for r in records:
        t, md = r["t"], r["mean_d"]
        mask  = t <= window_s
        if mask.sum() < 5: continue
        clipped.append((t[mask], md[mask] - md[mask][0]))

    if not clipped:
        print("  No recordings long enough for early drift — skipping."); return

    t_end = min(tc[-1] for tc, _ in clipped)
    t_com = np.linspace(0, t_end, 500)
    mat   = np.array([np.interp(t_com, tc, mdc) for tc, mdc in clipped])
    mu, sd = mat.mean(0), mat.std(0)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    for tc, mdc in clipped:
        axes[0].plot(tc, mdc, color="k", lw=0.6, alpha=0.2)
    axes[0].axhline(0, color="gray", lw=0.5, ls="--")
    axes[0].set_ylabel("Motion from onset (µm)")
    axes[0].set_title(f"Early drift — all recordings (n={len(clipped)})")
    _despine_top_right(axes[0])

    axes[1].fill_between(t_com, mu - sd, mu + sd, color="steelblue", alpha=0.25)
    axes[1].plot(t_com, mu, color="steelblue", lw=2)
    axes[1].axhline(0, color="gray", lw=0.5, ls="--")
    axes[1].set_ylabel("Motion from onset (µm)")
    axes[1].set_xlabel("Time from session start (s)")
    axes[1].set_title("Mean ± std across recordings")
    _despine_top_right(axes[1])

    fig.suptitle(f"Drift in first {window_s//60} minutes")
    fig.tight_layout()
    fig.savefig(out / "L2M_early_drift.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def _psth_matrices(records, psth_key, t_win_key):
    """Collect per-recording mean PSTH traces onto a common t_win grid."""
    rec_means, t_win_ref = [], None
    for r in records:
        psth = r.get(psth_key)
        t_win = r.get(t_win_key)
        if psth is None or t_win is None or len(psth) == 0:
            continue
        rec_means.append(psth.mean(axis=0))
        if t_win_ref is None:
            t_win_ref = t_win
    if not rec_means:
        return None, None
    # interpolate all onto the first recording's t_win grid
    mat = np.array([np.interp(t_win_ref, r.get(t_win_key, t_win_ref), rm)
                    if r.get(t_win_key) is not None else rm
                    for r, rm in zip(
            [r for r in records if r.get(psth_key) is not None
             and r.get(t_win_key) is not None
             and len(r.get(psth_key, [])) > 0],
            rec_means)])
    return mat, t_win_ref

def _plot_psth_figures(mat, t_win, n_trials_label, stim_label, color, out, prefix):
    grand = mat.mean(axis=0)
    sem = mat.std(axis=0) / np.sqrt(len(mat))

    # overlay
    fig, ax = plt.subplots(figsize=(5, 5))
    for trace in mat:
        ax.plot(t_win, trace, color="k", lw=0.7, alpha=0.2)
    ax.fill_between(t_win, grand - sem, grand + sem, color="k", alpha=0.2)
    ax.plot(t_win, grand, color="k", lw=2, label=f"mean ± SEM (n={len(mat)} recordings)")
    ax.axvline(0, color="k", lw=1.0, ls="--", label=stim_label)
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Time (s)");
    ax.set_ylabel("Motion from baseline (µm)")
    ax.set_title(f"{stim_label} — motion PSTH (interpolated)")
    ax.legend(frameon=False, fontsize=7)
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / f"{prefix}_overlay.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

def plot_whisker_psth(records, out):
    mat, t_win = _psth_matrices(records, "whisker_psth", "whisker_t_win")
    if mat is None:
        print("  No whisker PSTH data — skipping.");
        return
    _plot_psth_figures(mat, t_win, len(mat), "whisker stim", "steelblue",
                       out, "L2N_whisker_psth")

def plot_valve_psth(records, out):
    mat, t_win = _psth_matrices(records, "valve_psth", "valve_t_win")
    if mat is None:
        print("  No valve PSTH data — skipping.");
        return
    _plot_psth_figures(mat, t_win, len(mat), "valve open", "darkorange",
                       out, "L2O_valve_psth")


def plot_whisker_psth_old(records, out, pre: float = 1.0, post: float = 2.0, n_pts: int = 150):
    """Motion PSTH around whisker stimulation. Fig A: overlay. Fig B: mean ± SEM."""
    t_win     = np.linspace(-pre, post, n_pts)
    rec_means = []
    labels    = []
    for r in records:
        psth_data = r.get("whisker_psth")
        if psth_data is None or len(psth_data) == 0: continue
        psth, t_win = psth_data  # unpack tuple
        rec_means.append(psth.mean(axis=0))
        labels.append(f"{r['mouse']} | {r['session'][-10:]} | imec{r['probe']}")

    if not rec_means:
        print("  No whisker PSTH data — skipping."); return

    mat   = np.array(rec_means)
    grand = mat.mean(axis=0)
    sem   = mat.std(axis=0) / np.sqrt(len(mat))

    # Fig A: overlay
    fig, ax = plt.subplots(figsize=(8, 4))
    for trace in mat:
        ax.plot(t_win, trace, color="k", lw=0.7, alpha=0.2)
    ax.fill_between(t_win, grand - sem, grand + sem, color="k", alpha=0.2)
    ax.plot(t_win, grand, color="k", lw=2, label=f"mean ± SEM (n={len(mat)})")
    ax.axvline(0, color="firebrick", lw=1.0, ls="--", label="stim onset")
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Time from whisker stim (s)")
    ax.set_ylabel("Motion from baseline (µm)")
    ax.set_title("Whisker stim — motion PSTH (overlay)")
    ax.legend(frameon=False, fontsize=7)
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2N_whisker_psth_overlay.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    # Fig B: mean ± SEM with individual probe means as thin lines
    fig, ax = plt.subplots(figsize=(7, 4))
    for trace in mat:
        ax.plot(t_win, trace, color="0.75", lw=0.7, alpha=0.6)
    ax.fill_between(t_win, grand - sem, grand + sem, color="steelblue", alpha=0.25)
    ax.plot(t_win, grand, color="steelblue", lw=2.5,
            label=f"mean ± SEM (n={len(mat)} recordings)")
    ax.axvline(0, color="firebrick", lw=1.0, ls="--", label="stim onset")
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Time from whisker stim (s)")
    ax.set_ylabel("Motion from baseline (µm)")
    ax.set_title("Whisker stim — motion PSTH (average)")
    ax.legend(frameon=False, fontsize=7)
    _despine_top_right(ax)
    fig.tight_layout()
    fig.savefig(out / "L2N_whisker_psth_average.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_valve_psth_old(records, out, pre: float = 1.0, post: float = 2.0, n_pts: int = 150):
    """Motion PSTH around valve opening times. Fig A: overlay. Fig B: mean ± SEM."""
    t_win     = np.linspace(-pre, post, n_pts)
    rec_means = []
    for r in records:
        psth = r.get("valve_psth")
        if psth is None or len(psth) == 0: continue
        rec_means.append(psth.mean(axis=0))

    if not rec_means:
        print("  No valve PSTH data — skipping."); return

    mat   = np.array(rec_means)
    grand = mat.mean(axis=0)
    sem   = mat.std(axis=0) / np.sqrt(len(mat))

    for suffix, extra in [("overlay", True), ("average", False)]:
        fig, ax = plt.subplots(figsize=(8, 4))
        if extra:
            for trace in mat:
                ax.plot(t_win, trace, color="k", lw=0.7, alpha=0.2)
            ax.fill_between(t_win, grand - sem, grand + sem, color="k", alpha=0.2)
            ax.plot(t_win, grand, color="k", lw=2, label=f"mean ± SEM (n={len(mat)})")
        else:
            for trace in mat:
                ax.plot(t_win, trace, color="0.75", lw=0.7, alpha=0.6)
            ax.fill_between(t_win, grand - sem, grand + sem, color="darkorange", alpha=0.25)
            ax.plot(t_win, grand, color="darkorange", lw=2.5,
                    label=f"mean ± SEM (n={len(mat)} recordings)")
        ax.axvline(0, color="firebrick", lw=1.0, ls="--", label="valve open")
        ax.axhline(0, color="gray", lw=0.5, ls="--")
        ax.set_xlabel("Time from valve open (s)")
        ax.set_ylabel("Motion from baseline (µm)")
        ax.set_title(f"Valve — motion PSTH ({suffix})")
        ax.legend(frameon=False, fontsize=7)
        _despine_top_right(ax)
        fig.tight_layout()
        fig.savefig(out / f"L2O_valve_psth_{suffix}.png", dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
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


def report_extreme_recordings(records):
    df = _scalar_df(records)
    metrics = [
        ("max_abs_disp_um",   "Max |disp|"),
        ("p2p_um",            "P2P"),
        ("rms_um",            "RMS"),
        ("mean_abs_vel_um_s", "Mean |vel|"),
        ("n_fast_events_permissive", "Fast events"),
        ("tau",               "Settling τ"),
        ("tau_r2",            "Tau R²"),
    ]
    print("\n── Extreme recordings (top/bottom 10% per metric) ──────────────────")
    flagged = defaultdict(list)

    for metric, label in metrics:
        if metric not in df.columns: continue
        col = df[metric].dropna()
        if len(col) < 5: continue
        lo, hi = col.quantile(0.10), col.quantile(0.90)
        for _, row in df.iterrows():
            val = row.get(metric)
            if pd.isna(val): continue
            key = (row["mouse"], row["session"], int(row["probe"]))
            if val <= lo:
                flagged[key].append(f"{label}={val:.3g} [bottom 10%, thr≤{lo:.3g}]")
            elif val >= hi:
                flagged[key].append(f"{label}={val:.3g} [top 10%, thr≥{hi:.3g}]")

    if not flagged:
        print("  None found."); return
    for (mouse, session, probe), reasons in sorted(flagged.items()):
        print(f"\n  {mouse} | {session} | imec{probe}")
        for r in reasons:
            print(f"    · {r}")
    print()

def _load_and_compute(args):
    mouse, session, probe, dredge_path, reward_group, insertion_depths, target_areas, day_of_recording = args
    global INSERTION_DEPTHS, TARGET_AREAS, DAY_OF_RECORDING
    INSERTION_DEPTHS  = insertion_depths
    TARGET_AREAS      = target_areas
    DAY_OF_RECORDING  = day_of_recording
    try:
        motion = load_motion_info(Path(dredge_path))["motion"]
    except FileNotFoundError as e:
        return "SKIP", mouse, session, probe, None, str(e)
    if motion is None:
        return "SKIP", mouse, session, probe, None, "no motion object"
    try:
        metrics = compute_metrics(mouse, reward_group, session, probe, motion,
                                  dredge_path=Path(dredge_path))
    except Exception as e:
        return "ERR", mouse, session, probe, None, str(e)
    return "OK", mouse, session, probe, metrics, None

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Probe drift motion summary.")
    parser.add_argument("--mode", choices=["single", "dataset", "all"], default="dataset") #["single", "dataset", "all"]
    parser.add_argument("--mouse",   default=None)
    parser.add_argument("--session", default=None)
    parser.add_argument("--workers", type=int, default=N_WORKERS)
    args = parser.parse_args()

    global INSERTION_DEPTHS, TARGET_AREAS, DAY_OF_RECORDING
    INSERTION_DEPTHS = _load_insertion_depths()
    REWARD_GROUPS    = _load_reward_group()
    TARGET_AREAS     = _load_target_area()
    DAY_OF_RECORDING = _load_day_of_recording()

    L1_OUT.mkdir(parents=True, exist_ok=True)
    L2_OUT.mkdir(parents=True, exist_ok=True)

    all_recs = list(discover_recordings(DATA_ROOT, args.mouse, args.session))
    if not all_recs:
        print("No recordings found."); return

    records:  list[dict]      = []
    by_mouse: dict[str, list] = defaultdict(list)

    #print(f"Loading {len(all_recs)} recording(s)...")
    #for mouse, session, probe, dredge_path in all_recs:
    #    reward_group = REWARD_GROUPS.get(mouse, "unknown")
    #    try:
    #        motion = load_motion_info(dredge_path)["motion"]
    #    except FileNotFoundError as e:
    #        print(f"  SKIP {mouse}/{session}/imec{probe}: {e}"); continue
    #    if motion is None:
    #        print(f"  SKIP {mouse}/{session}/imec{probe}: no motion object"); continue
    #    by_mouse[mouse].append((session, probe, motion))
    #    records.append(compute_metrics(mouse, reward_group, session, probe, motion,
    #                                   dredge_path=dredge_path))
    #    print(f"  loaded {mouse} | {session} | imec{probe}")

    print(f"Loading {len(all_recs)} recording(s)...")
    worker_args = [
        (mouse, session, probe, str(dredge_path),
         REWARD_GROUPS.get(mouse, "unknown"),
         INSERTION_DEPTHS, TARGET_AREAS, DAY_OF_RECORDING)
        for mouse, session, probe, dredge_path in all_recs
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_load_and_compute, a): a for a in worker_args}
        for fut in as_completed(futs):
            status, mouse, session, probe, metrics, msg = fut.result()
            if status != "OK":
                print(f"  {status} {mouse}/{session}/imec{probe}: {msg}")
                continue
            records.append(metrics)
            print(f"  loaded {mouse} | {session} | imec{probe}")

    records.sort(key=lambda r: (r["mouse"], r["session"], r["probe"]))

    # Reload motion objects for by_mouse — only needed for L1
    if args.mode in ("single", "all"):
        print("Reloading motion objects for L1 figures...")
        for mouse, session, probe, dredge_path in all_recs:
            try:
                motion = load_motion_info(dredge_path)["motion"]
            except FileNotFoundError:
                continue
            if motion is None:
                continue
            by_mouse[mouse].append((session, probe, motion))

    if not records:
        print("No valid motion objects found."); return

    if args.mode in ("single", "all"):
        print(f"\nGenerating L1 per-recording figures ({args.workers} workers)...")
        worker_args = [
            (mouse, session, probe, str(dredge_path), str(L1_OUT), INSERTION_DEPTHS, TARGET_AREAS)
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
        plot_fast_event_summary(records, L2_OUT)
        plot_drift_over_session(records, L2_OUT)
        plot_exponential_decay(records, L2_OUT)
        plot_autocorrelation(records, L2_OUT)
        plot_depth_gradient(records, L2_OUT)
        plot_depth_gradient_absolute(records, L2_OUT)
        plot_depth_gradient_barh(records, L2_OUT)
        plot_insertion_depth_scatter(records, L2_OUT)
        plot_motion_by_target(records, L2_OUT)
        plot_learning_vs_expert(records, L2_OUT)
        plot_early_drift(records, L2_OUT)
        plot_whisker_psth(records, L2_OUT)
        plot_valve_psth(records, L2_OUT)
        report_extreme_recordings(records)

    print(f"\nDone.\n  L1 : {L1_OUT}\n  L2 : {L2_OUT}")


if __name__ == "__main__":
    main()
