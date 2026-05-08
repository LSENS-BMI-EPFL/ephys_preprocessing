"""
Probe drift motion summary across Neuropixels recordings.

Path conventions
----------------
Per-recording figures : L1_OUT/<mouse>_<session>_imecN.png
Per-mouse figures     : L1_OUT/<mouse>_all_probes.png
Dataset figures/table : L2_OUT
"""

from __future__ import annotations
from collections import defaultdict
import json
from pathlib import Path
import warnings

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.signal import find_peaks

from spikeinterface.core.motion import Motion


def _despine_top_right(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _despine_top_right_bottom(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)

def load_motion_info(folder):
    folder = Path(folder)
    motion_info = {}
    with open(folder / "parameters.json") as f:
        motion_info["parameters"] = json.load(f)
    with open(folder / "run_times.json") as f:
        motion_info["run_times"] = json.load(f)
    for name in ("peaks", "peak_locations"):
        p = folder / f"{name}.npy"
        motion_info[name] = np.load(p) if p.exists() else None
    if (folder / "motion").is_dir():
        motion = Motion.load(folder / "motion")
    else:
        required = ["spatial_bins.npy", "temporal_bins.npy", "motion.npy"]
        if all((folder / f).is_file() for f in required):
            warnings.warn("Trying to load Motion from the legacy format")
            motion = Motion(
                displacement=[np.load(folder / "motion.npy")],
                temporal_bins_s=[np.load(folder / "temporal_bins.npy")],
                spatial_bins_um=np.load(folder / "spatial_bins.npy"),
            )
        else:
            warnings.warn("No `motion` object in folder. `motion` set to None.")
            motion = None
    motion_info["motion"] = motion
    return motion_info

# ── Config ─────────────────────────────────────────────────────────────────
DATA_ROOT = Path(r"M:\analysis\Axel_Bisi\data")
COMBINED  = Path(r"M:\analysis\Axel_Bisi\combined_results\motion_estimation")
L1_OUT    = COMBINED / "single"
L2_OUT    = COMBINED

# Fast-event detection: (permissive_thr, conservative_thr) in µm/s
THRESHOLDS = (5.0, 10.0)
THR_COLORS = ("darkorange", "firebrick")
THR_LABELS = [f"permissive (>{THRESHOLDS[0]} µm/s)", f"conservative (>{THRESHOLDS[1]} µm/s)"]

CMAP_DISP      = "PiYG"
FIG_DPI        = 150
INSERTION_INFO = Path(r"M:\share_internal\Axel_Bisi_Share\dataset_info\joint_probe_insertion_info.xlsx")

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8})


def _load_insertion_depths() -> dict[tuple[str, str], float]:
    """Return {(mouse_name, probe_id): depth_um}. Skips if file missing."""
    if not INSERTION_INFO.exists():
        return {}
    df = pd.read_excel(INSERTION_INFO)
    return {(str(row.mouse_name), str(row.probe_id)): float(row.depth)
            for row in df.itertuples(index=False)}


INSERTION_DEPTHS: _load_insertion_depths()


# ── Path helpers ────────────────────────────────────────────────────────────

def discover_recordings(root: Path):
    """Yield (mouse, session, probe_str, dredge_fast_path) for every dredge_fast folder."""
    for p in sorted(root.glob("*/*/Ephys/catgt_*/*_imec*/dredge_fast/motion")):
        parts = p.parts
        base    = len(root.parts)
        mouse   = parts[base]
        session = parts[base + 1]
        probe   = int(parts[-3][-1:])
        yield mouse, session, probe, p


def l1_dir() -> Path:
    L1_OUT.mkdir(parents=True, exist_ok=True)
    return L1_OUT


# ── Signal helpers ──────────────────────────────────────────────────────────

def mean_disp(motion: Motion, seg: int = 0) -> np.ndarray:
    return np.mean(motion.displacement[seg], axis=1)


def vel(motion: Motion, seg: int = 0) -> np.ndarray:
    t = motion.temporal_bins_s[seg]
    return np.gradient(mean_disp(motion, seg), t)


def detect_peaks_dual(t: np.ndarray, velocity: np.ndarray):
    """Return (peaks_permissive_idx, peaks_conservative_idx)."""
    p_per, _ = find_peaks(np.abs(velocity), height=THRESHOLDS[0])
    p_con, _ = find_peaks(np.abs(velocity), height=THRESHOLDS[1])
    return p_per, p_con


def _filter_by_depth(disp: np.ndarray, dep: np.ndarray, mouse: str, probe: str):
    """Mask spatial bins to those shallower than insertion depth. Returns (disp, dep)."""
    depth_um = INSERTION_DEPTHS.get((mouse, str(probe)))
    if depth_um is not None:
        mask = dep < depth_um
        return disp[:, mask], dep[mask]
    return disp, dep


# ── Level 1-R : per-recording figure ────────────────────────────────────────

def plot_recording(mouse: str, session: str, probe: str, motion: Motion, out_dir: Path):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    seg  = 0
    disp = motion.displacement[seg]
    t    = motion.temporal_bins_s[seg]
    dep  = motion.spatial_bins_um
    disp, dep = _filter_by_depth(disp, dep, mouse, probe)
    md   = np.mean(disp, axis=1)
    v    = np.gradient(md, t)
    p_per, p_con = detect_peaks_dual(t, v)

    depth_um = INSERTION_DEPTHS.get((mouse, str(probe)))
    title = f"{mouse}  |  {session}  |  imec{probe}"
    if depth_um is not None:
        title += f" | insertion depth {depth_um:.0f} µm"

    fig = plt.figure(figsize=(13, 11), dpi=500)
    gs  = gridspec.GridSpec(4, 1, hspace=0.06, height_ratios=[2.2, 1.3, 1.3, 1.5],
                            top=0.94, bottom=0.07, left=0.08, right=0.93)
    ax_h = fig.add_subplot(gs[0])
    ax_d = fig.add_subplot(gs[1], sharex=ax_h)
    ax_v = fig.add_subplot(gs[2], sharex=ax_h)
    ax_e = fig.add_subplot(gs[3], sharex=ax_h)

    # ── P1 heatmap — horizontal colorbar via divider to preserve x-alignment ─
    clim = np.max(np.abs(disp)) * 1.05
    im = ax_h.imshow(disp.T, aspect="auto", origin="lower", cmap=CMAP_DISP,
                     extent=(t[0], t[-1], dep[0], dep[-1]),
                     vmin=-clim, vmax=clim, interpolation="nearest")
    divider = make_axes_locatable(ax_h)
    cax = divider.append_axes("bottom", size="6%", pad=0.05)
    cax.set_axis_off()

    # create smaller centered axes inside it
    cbax = inset_axes(
        cax,
        width="25%",  # colorbar width relative to cax
        height="80%",
        loc="center",
    )

    cb = fig.colorbar(im, cax=cbax, orientation="horizontal", aspect=30)

    cb.set_label("Motion (µm)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ## shrink the colorbar axes manually
    #pos = cax.get_position() # shrink colorbar manually
    #new_width = pos.width * 0.5  # 50% width
    #new_x = pos.x0 + pos.width * 0.25  # center it
    #cax.set_position([new_x, pos.y0, new_width, pos.height])
    #cb = fig.colorbar(im, ax=cax, orientation="horizontal", shrink=0.5)
    #cb.set_label("Motion (µm)", fontsize=7)
    #cb.ax.tick_params(labelsize=6)
    ax_h.set_ylabel("Depth (µm)")
    ax_h.set_title(title)
    plt.setp(ax_h.get_xticklabels(), visible=False)


    # ── P2 traces + mean — coloured by depth, with legend ───────────────────
    n_d    = disp.shape[1]
    cmap_d = plt.colormaps["viridis"].resampled(max(n_d, 1))
    # Show at most ~8 legend entries to avoid clutter; pick evenly spaced indices
    max_legend = 8
    legend_idx = np.round(np.linspace(0, n_d - 1, min(n_d, max_legend))).astype(int)
    legend_set = set(legend_idx)
    for i in range(n_d):
        lbl = f"{dep[i]:.0f} µm" if (i in legend_set) else None
        ax_d.plot(t, disp[:, i], lw=0.8, alpha=0.45, color=cmap_d(i), label=lbl)
    ax_d.plot(t, md, color="k", lw=1.4, label="mean")
    ax_d.axhline(0, color="gray", lw=0.5, ls="--")
    ax_d.set_ylabel("Motion (µm)")
    ax_d.legend(loc="upper right", fontsize=6, handlelength=1,
                ncol=2, frameon=False)
    plt.setp(ax_d.get_xticklabels(), visible=False)

    # ── P3 velocity ──────────────────────────────────────────────────────────
    ax_v.plot(t, v, color="steelblue", lw=0.9)
    ax_v.axhline(0, color="gray", lw=0.5, ls="--")
    ax_v.set_ylabel("Velocity (µm/s)")
    plt.setp(ax_v.get_xticklabels(), visible=False)

    # ── P4 fast events ────────────────────────────────────────────────────────
    ax_e.plot(t, v, color="steelblue", lw=0.85, alpha=0.6, label="velocity")
    ax_e2 = ax_e.twinx()
    ax_e2.plot(t, md, color="0.25", lw=1.1, alpha=0.8, label="motion")
    ax_e2.set_ylabel("Motion (µm)", color="0.25")

    for peaks_idx, color, label in zip((p_per, p_con), THR_COLORS, THR_LABELS):
        ax_e.scatter(t[peaks_idx], v[peaks_idx], color=color, s=35, zorder=5,
                     label=f"{label}  n={len(peaks_idx)}", linewidths=0)
        for pi in peaks_idx:
            ax_e.axvline(t[pi], color=color, lw=0.7, alpha=0.5)

    for thr, color in zip(THRESHOLDS, THR_COLORS):
        ax_e.axhline( thr, color=color, lw=0.7, ls=":", alpha=0.8)
        ax_e.axhline(-thr, color=color, lw=0.7, ls=":", alpha=0.8)

    ax_e.set_ylabel("Velocity (µm/s)")
    ax_e.set_xlabel("Time (s)")
    handles_v, labels_v = ax_e.get_legend_handles_labels()
    handles_m, labels_m = ax_e2.get_legend_handles_labels()

    ax_e.legend(handles_v + handles_m, labels_v + labels_m, fontsize=6.5,
                loc="lower center", bbox_to_anchor=(0.5, -0.35),
                ncol=4, handlelength=1, frameon=False)

    for idx, _ax in enumerate((ax_d, ax_v, ax_e, ax_e2)):
        _despine_top_right_bottom(_ax)

    # Align x and y labels
    fig.tight_layout()
    fig.align_labels()
    fig.subplots_adjust(hspace=0.7)
    rec_label = f"{mouse}_{session}_imec{probe}"
    fig.savefig(out_dir / f"{rec_label}.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── Level 1-M : per-mouse figure ────────────────────────────────────────────

def plot_mouse(mouse: str, recs: list[tuple[str, str, Motion]], out_dir: Path):
    """recs : list of (session, probe, motion)"""
    n = len(recs)
    colors = plt.colormaps["tab10"].resampled(max(n, 1))

    fig = plt.figure(figsize=(13, 10))
    gs  = gridspec.GridSpec(3, 1, hspace=0.35, height_ratios=[2, 2, 2.5],
                            top=0.93, bottom=0.07)
    ax_d = fig.add_subplot(gs[0])
    ax_v = fig.add_subplot(gs[1], sharex=ax_d)
    ax_c = fig.add_subplot(gs[2])

    series, plot_labels = [], []
    for i, (session, probe, motion) in enumerate(recs):
        t    = motion.temporal_bins_s[0]
        disp, _ = _filter_by_depth(motion.displacement[0], motion.spatial_bins_um, mouse, probe)
        md   = np.mean(disp, axis=1)
        v    = np.gradient(md, t)
        lbl = f"imec{probe} | {session}"
        plot_labels.append(lbl)
        series.append((t, md))
        c = colors(i)
        ax_d.plot(t, md, color=c, lw=1.2, label=lbl)
        ax_v.plot(t, v,  color=c, lw=0.9, alpha=0.85)

    ax_d.axhline(0, color="gray", lw=0.5, ls="--")
    ax_d.set_ylabel("Mean motion (µm)")
    ax_d.set_title(f"{mouse} — all probes")
    ax_d.legend(fontsize=6.5, loc="upper right", ncol=min(3, n))
    plt.setp(ax_d.get_xticklabels(), visible=False)

    ax_v.axhline(0, color="gray", lw=0.5, ls="--")
    ax_v.set_ylabel("Velocity (µm/s)")
    ax_v.set_xlabel("Time (s)")


    # Pearson correlation matrix on common time grid
    t_end = min(s[0][-1] for s in series)
    t_com = np.linspace(0, t_end, 500)
    mat   = np.array([np.interp(t_com, s[0], s[1]) for s in series])
    corr  = np.corrcoef(mat)

    im = ax_c.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax_c.set_xticks(range(n)); ax_c.set_xticklabels(plot_labels, rotation=30, ha="right", fontsize=7)
    ax_c.set_yticks(range(n)); ax_c.set_yticklabels(plot_labels, fontsize=7)
    for i in range(n):
        for j in range(n):
            ax_c.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=7,
                      color="w" if abs(corr[i,j]) > 0.65 else "k")
    plt.colorbar(im, ax=ax_c, label="Pearson r", pad=0.01, fraction=0.025)
    ax_c.set_title("Inter-probe displacement correlation", fontsize=8)

    fig.savefig(out_dir / f"{mouse}_all_probes.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── Metrics ──────────────────────────────────────────────────────────────────

def compute_metrics(mouse: str, session: str, probe: str, motion: Motion) -> dict:
    t    = motion.temporal_bins_s[0]
    disp, _ = _filter_by_depth(motion.displacement[0], motion.spatial_bins_um, mouse, probe)
    md   = np.mean(disp, axis=1)
    v    = np.gradient(md, t)
    p_per, p_con = detect_peaks_dual(t, v)
    return dict(
        mouse=mouse, session=session, probe=probe,
        duration_s=float(t[-1]),
        max_abs_disp_um=float(np.max(np.abs(md))),
        p2p_um=float(md.max() - md.min()),
        rms_um=float(np.sqrt(np.mean(md ** 2))),
        mean_abs_vel_um_s=float(np.mean(np.abs(v))),
        n_fast_events_conservative=len(p_con),
        n_fast_events_permissive=len(p_per),
        t=t, mean_d=md,
    )


# ── Level 2 figures ──────────────────────────────────────────────────────────

def _mouse_colors(records: list[dict]) -> dict:
    mice = sorted(set(r["mouse"] for r in records))
    cmap = plt.colormaps["tab10"].resampled(max(len(mice), 1))
    return {m: cmap(i) for i, m in enumerate(mice)}


def _interp_common(records: list[dict], n: int = 500):
    t_end = min(r["t"][-1] for r in records)
    t_com = np.linspace(0, t_end, n)
    mat   = np.array([np.interp(t_com, r["t"], r["mean_d"]) for r in records])
    return t_com, mat


def plot_grand_overlay(records: list[dict], out: Path):
    # stagger traces vertically by a fixed offset derived from overall amplitude
    all_amp = np.concatenate([r["mean_d"] for r in records])
    spacing = np.percentile(np.abs(all_amp), 90) * 2.5

    fig, ax = plt.subplots(figsize=(12, max(5, len(records) * 0.6)))
    yticks, ylabels = [], []
    for i, r in enumerate(records):
        offset = i * spacing
        ax.plot(r["t"], r["mean_d"] + offset, color="k", lw=0.8, alpha=0.75)
        yticks.append(offset)
        ylabels.append(f"{r['mouse']} | {r['session'][-10:]} | imec{r['probe']}")

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean motion (µm)  (+ offset per recording)")
    ax.set_title("Grand overlay — all recordings")
    fig.tight_layout()
    fig.savefig(out / "L2A_grand_overlay.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_magnitude_distributions(records: list[dict], out: Path):
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("t", "mean_d")}
                       for r in records])
    metrics = ["max_abs_disp_um", "p2p_um", "rms_um", "mean_abs_vel_um_s"]
    ylabels = ["Max |disp| (µm)", "Peak-to-peak (µm)", "RMS (µm)", "Mean |vel| (µm/s)"]

    fig, axes = plt.subplots(4, 1, figsize=(8, 14), sharex=False)
    for ax, metric, ylabel in zip(axes, metrics, ylabels):
        sns.violinplot(data=df, x="mouse", y=metric, ax=ax, inner=None,
                       color="0.75", linewidth=0.8)
        sns.stripplot(data=df, x="mouse", y=metric, ax=ax, color="k",
                      size=4, jitter=True, linewidth=0)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=30)
    axes[-1].set_xlabel("Mouse")
    fig.suptitle("Motion magnitude distributions")
    fig.tight_layout()
    fig.savefig(out / "L2B_magnitude_distributions.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_cross_recording_correlation(records: list[dict], out: Path):
    t_com, mat = _interp_common(records)
    n    = len(records)
    mice = [r["mouse"] for r in records]
    mc   = _mouse_colors(records)

    # collect all pair data
    all_x, all_y = [], []
    same_x, same_y, same_c = [], [], []
    diff_x, diff_y = [], []

    for i in range(n):
        for j in range(i + 1, n):
            x, y = mat[i], mat[j]
            all_x.extend(x); all_y.extend(y)
            if mice[i] == mice[j]:
                same_x.extend(x); same_y.extend(y)
                same_c.extend([mc[mice[i]]] * len(x))
            else:
                diff_x.extend(x); diff_y.extend(y)

    scatter_kw = dict(s=1, alpha=0.08, linewidths=0, rasterized=True)

    # Fig 1: all pairs in black
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(all_x, all_y, color="k", **scatter_kw)
    ax.set_xlabel("Motion recording A (µm)")
    ax.set_ylabel("Motion recording B (µm)")
    ax.set_title("Pairwise motion — all pairs")
    ax.set_aspect("equal", "datalim")
    fig.tight_layout()
    fig.savefig(out / "L2D_pairwise_scatter_all.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    # Fig 2: same-mouse colored, cross-mouse gray
    fig, ax = plt.subplots(figsize=(5, 5))
    if diff_x:
        ax.scatter(diff_x, diff_y, color="0.75", **scatter_kw)
    if same_x:
        ax.scatter(same_x, same_y, c=same_c, **scatter_kw)
    ax.set_xlabel("Motion recording A (µm)")
    ax.set_ylabel("Motion recording B (µm)")
    ax.set_title("Pairwise motion — same mouse colored")
    ax.set_aspect("equal", "datalim")
    fig.tight_layout()
    fig.savefig(out / "L2D_pairwise_scatter_bymouse.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_fast_event_summary(records: list[dict], out: Path):
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("t", "mean_d")}
                       for r in records])
    df = df.sort_values(["mouse", "session", "probe"]).reset_index(drop=True)
    lbls = [f"{r.mouse} | {r.session[-15:]} | imec{r.probe}" for r in df.itertuples()]
    y    = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.45)))
    ax.barh(y, df["n_fast_events_permissive"],   color=THR_COLORS[1], alpha=0.55, label=THR_LABELS[1])
    ax.barh(y, df["n_fast_events_conservative"], color=THR_COLORS[0], alpha=0.9,  label=THR_LABELS[0])
    ax.set_yticks(y); ax.set_yticklabels(lbls, fontsize=7)
    ax.set_xlabel("Number of fast events")
    ax.set_title("Fast-event count per recording")
    # Make legend in bottom of the plot to avoid covering bars; use smaller font to fit long labels
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "L2E_fast_event_summary.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_drift_over_session(records: list[dict], out: Path):
    by_session: dict[tuple, list] = defaultdict(list)
    for r in records:
        by_session[(r["mouse"], r["session"])].append(r)
    for k in by_session:
        by_session[k].sort(key=lambda r: r["probe"])

    fig, ax = plt.subplots(figsize=(7, 5))
    for (mouse, session), recs in by_session.items():
        xs = [i + 1 for i in range(len(recs))]
        ys = [r["mean_abs_vel_um_s"] for r in recs]
        ax.plot(xs, ys, color="0.5", lw=0.9, alpha=0.6)
        ax.scatter(xs, ys, color="k", s=25, zorder=3, linewidths=0)

    ax.set_xlabel("Probe index within session")
    ax.set_ylabel("Mean |velocity| (µm/s)")
    ax.set_title("Probe drift across session")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.tight_layout()
    fig.savefig(out / "L2C_drift_over_session.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def save_summary_table(records: list[dict], out: Path):
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("t", "mean_d")}
                       for r in records])
    df = df.sort_values(["mouse", "session", "probe"]).reset_index(drop=True)
    df.to_csv(out / "motion_summary.csv", index=False)

    # rendered PNG table
    numeric_cols = df.select_dtypes(include="number").columns
    df[numeric_cols] = df[numeric_cols].round(3)

    fig, ax = plt.subplots(figsize=(16, max(2, len(df) * 0.38 + 0.8)))
    ax.axis("off")
    tbl = ax.table(cellText=df.values.tolist(), colLabels=df.columns.tolist(),
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6.5)
    tbl.auto_set_column_width(range(len(df.columns)))
    fig.tight_layout()
    fig.savefig(out / "L2F_summary_table.png", dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    global INSERTION_DEPTHS
    INSERTION_DEPTHS = _load_insertion_depths()
    L2_OUT.mkdir(parents=True, exist_ok=True)
    all_recs = list(discover_recordings(DATA_ROOT))
    if not all_recs:
        print("No recordings found under", DATA_ROOT)
        return

    records:       list[dict]              = []
    by_mouse: dict[str, list]              = defaultdict(list)

    for mouse, session, probe, dredge_path in all_recs:
        print(f"  {mouse} | {session} | imec{probe}")
        try:
            motion = load_motion_info(dredge_path)["motion"]
        except FileNotFoundError:
            print(f"No motion estimation found")

        if motion is None:
            print(f"    ↳ skipped (no motion object)")
            continue

        #plot_recording(mouse, session, probe, motion, l1_dir())

        by_mouse[mouse].append((session, probe, motion))
        records.append(compute_metrics(mouse, session, probe, motion))

    for mouse, recs in by_mouse.items():
        plot_mouse(mouse, recs, l1_dir())

    plot_grand_overlay(records, L2_OUT)
    plot_magnitude_distributions(records, L2_OUT)
    plot_cross_recording_correlation(records, L2_OUT)
    plot_fast_event_summary(records, L2_OUT)
    plot_drift_over_session(records, L2_OUT)
    save_summary_table(records, L2_OUT)

    print("\nDone.")
    print(f"  L1 (single)  : {L1_OUT}")
    print(f"  L2 (combined): {L2_OUT}")


if __name__ == "__main__":
    main()
