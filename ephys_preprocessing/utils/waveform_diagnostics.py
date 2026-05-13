"""
waveform_diagnostics.py
---------------------
Single-figure diagnostic summary for C_Waves outputs.

Expected files (set --data_dir, or place script alongside outputs):
  - mean_waveforms.npy          (nClusters, nChannels, nSamples)  µV
  - median_peak_waveforms.npy   (nClusters, nSamples)             µV
  - cluster_snr.npy             (nClusters, 2)  -> [snr, nSpikes_in_snr]
  - waveform_metrics.csv  OR  mean_waveform.csv
        columns: cluster_id, peak_channel, duration, halfwidth,
                 pt_ratio, repolarization_slope, recovery_slope,
                 trough_idx, peak_idx

Usage:
  python waveform_diagnostics.py --data_dir /path/to/cwaves/output
  python waveform_diagnostics.py --data_dir . --out diagnostics.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── dark theme ─────────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "figure.facecolor": "#0d1117",
    "axes.facecolor":   "#161b22",
    "axes.edgecolor":   "#30363d",
    "axes.labelcolor":  "#c9d1d9",
    "axes.titlecolor":  "#e6edf3",
    "xtick.color":      "#8b949e",
    "ytick.color":      "#8b949e",
    "text.color":       "#c9d1d9",
    "grid.color":       "#21262d",
    "grid.linewidth":   0.6,
    "legend.facecolor": "#161b22",
    "legend.edgecolor": "#30363d",
    "font.size":        9,
})

BLUE   = "#58a6ff"
PURPLE = "#bc8cff"
GOLD   = "#e3b341"
GREEN  = "#3fb950"
RED    = "#f85149"
DIM    = "#8b949e"

SNR_THRESHOLDS = {"poor": 2.0, "ok": 4.0}
SNR_VMIN, SNR_VMAX = 0, 10
CMAP = "viridis"


# ── loading ────────────────────────────────────────────────────────────────────

def load_data(data_dir: Path) -> dict:
    d = {}

    for key, fname in [
        ("snr",     "cluster_snr.npy"),
        ("mean_wf", "mean_waveforms.npy"),
        ("med_wf",  "median_peak_waveforms.npy"),
    ]:
        p = data_dir / fname
        if not p.exists():
            sys.exit(f"[ERROR] Required file not found: {p}")
        d[key] = np.load(p)
        print(f"[OK] Loaded {fname}  shape={d[key].shape}")

    d["n_clusters"], d["n_channels"], d["n_samples"] = d["mean_wf"].shape

    # Accept either filename
    for csv_name in ("waveform_metrics.csv", "mean_waveform.csv"):
        csv_path = data_dir / csv_name
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            df["cluster_id"] = df["cluster_id"].astype(int)
            d["metrics"] = df
            print(f"[OK] Loaded {csv_name}  shape={df.shape}  cols={list(df.columns)}")
            break
    else:
        print("[WARN] No waveform metrics CSV found — metric panels will be skipped.")
        d["metrics"] = None

    return d


def build_masks(d: dict):
    snr  = d["snr"][:, 0].astype(float)
    nspk = d["snr"][:, 1].astype(float)
    valid = snr > 0   # SNR == -1 marks empty clusters
    return snr, nspk, valid


def align_metrics(metrics_df: pd.DataFrame, valid: np.ndarray) -> pd.DataFrame:
    """Return metrics rows for valid clusters only."""
    valid_ids = set(np.where(valid)[0])
    return metrics_df[metrics_df["cluster_id"].isin(valid_ids)].copy()


# ── colour / annotation helpers ────────────────────────────────────────────────

def snr_colors(snr_vals):
    norm = Normalize(vmin=SNR_VMIN, vmax=SNR_VMAX)
    return plt.get_cmap(CMAP)(norm(snr_vals))


def add_snr_colorbar(ax, fig, label="SNR"):
    sm = ScalarMappable(cmap=CMAP, norm=Normalize(SNR_VMIN, SNR_VMAX))
    sm.set_array([])
    return fig.colorbar(sm, ax=ax, label=label, pad=0.02, fraction=0.046)


def _threshold_lines(ax, thresholds=SNR_THRESHOLDS, orient="v"):
    fn = ax.axvline if orient == "v" else ax.axhline
    fn(thresholds["poor"], color=RED,  lw=1.1, ls="--", alpha=0.8,
       label=f'poor ({thresholds["poor"]})')
    fn(thresholds["ok"],   color=GOLD, lw=1.1, ls="--", alpha=0.8,
       label=f'ok ({thresholds["ok"]})')


def _note(ax, txt, loc="ur"):
    xp, ha = (0.97, "right") if "r" in loc else (0.03, "left")
    yp, va = (0.97, "top")   if "u" in loc else (0.03, "bottom")
    ax.text(xp, yp, txt, transform=ax.transAxes, ha=ha, va=va,
            fontsize=7.5, color=DIM,
            bbox=dict(boxstyle="round,pad=0.3", fc="#0d1117", alpha=0.7))


def _no_data(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes,
            ha="center", va="center", color=DIM, fontsize=8)


# ── panels ─────────────────────────────────────────────────────────────────────

def panel_snr_hist(ax, snr, valid):
    snr_v = snr[valid]
    ax.hist(snr_v, bins=60, color=BLUE, alpha=0.85, edgecolor="none")
    _threshold_lines(ax, orient="v")
    med = np.median(snr_v)
    ax.axvline(med, color=GREEN, lw=1.4, label=f"median={med:.2f}")
    ax.set_xlabel("SNR")
    ax.set_ylabel("# clusters")
    ax.set_title("SNR Distribution")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y")
    poor = (snr_v < SNR_THRESHOLDS["poor"]).sum()
    _note(ax, f"n={valid.sum()} valid  |  {poor} poor", "ur")


def panel_nspikes_hist(ax, nspk, valid):
    ax.hist(nspk[valid], bins=50, color=PURPLE, alpha=0.85, edgecolor="none")
    med = np.median(nspk[valid])
    ax.axvline(med, color=GREEN, lw=1.4, label=f"median={med:.0f}")
    ax.set_xlabel("nSpikes used for SNR")
    ax.set_ylabel("# clusters")
    ax.set_title("Spike Count Distribution")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y")


def panel_snr_vs_nspk(ax, snr, nspk, valid):
    sc = ax.scatter(nspk[valid], snr[valid], c=snr[valid],
                    cmap=CMAP, vmin=SNR_VMIN, vmax=SNR_VMAX,
                    s=8, alpha=0.7, linewidths=0)
    plt.colorbar(sc, ax=ax, label="SNR", pad=0.02)
    _threshold_lines(ax, orient="h")
    ax.set_xlabel("nSpikes used for SNR")
    ax.set_ylabel("SNR")
    ax.set_title("SNR vs Spike Count")
    ax.legend(fontsize=7)
    ax.grid(True)
    ax.set_xscale("log")


def panel_peak_waveforms(ax, fig, med_wf, snr, valid, metrics_df, max_traces=20):
    """
    Overlay median peak waveforms coloured by SNR.
    x-axis is aligned to the median trough_idx so trough sits at sample 0.
    Vertical dashed lines mark the median trough and peak positions.
    """
    idx = np.where(valid)[0]
    if len(idx) > max_traces:
        # Sample evenly across SNR range so we see the full quality spectrum
        order = idx[np.argsort(snr[idx])]
        step  = max(1, len(order) // max_traces)
        idx   = order[::step][:max_traces]

    n_samples = med_wf.shape[1]

    # Reference point: median trough_idx from metrics (fallback: grand-mean trough)
    t_ref = 0
    med_peak_offset = None
    if metrics_df is not None and "trough_idx" in metrics_df.columns:
        df_v = align_metrics(metrics_df, valid)
        t_ref = int(df_v["trough_idx"].median())
        if "peak_idx" in df_v.columns:
            med_peak_offset = int(df_v["peak_idx"].median()) - t_ref
    else:
        grand = np.nanmean(med_wf[valid], axis=0)
        t_ref = int(np.argmin(grand))

    t = np.arange(n_samples) - t_ref

    # Draw traces
    colors = snr_colors(snr[idx])
    for i, c in zip(idx, colors):
        ax.plot(t, med_wf[i], color="white", lw=0.45, alpha=0.35)

    # Grand mean on top
    grand = np.nanmean(med_wf[valid], axis=0)
    ax.plot(t, grand, color="white", lw=1.6, label="grand mean", zorder=5)

    # Trough and peak reference lines
    ax.axvline(0, color=RED, lw=1.0, ls=":", alpha=0.9, label="median trough (t=0)")
    if med_peak_offset is not None:
        ax.axvline(med_peak_offset, color=GOLD, lw=1.0, ls=":", alpha=0.9,
                   label=f"median peak (+{med_peak_offset} samp)")

    ax.set_xlabel("Samples relative to trough")
    ax.set_ylabel("Amplitude (µV)")
    ax.set_title(f"Median Peak Waveforms  (n={len(idx)} shown, coloured by SNR)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True)
    add_snr_colorbar(ax, fig)


def panel_duration_halfwidth(ax, fig, metrics_df, snr, valid):
    """Scatter of duration vs halfwidth — separates narrow (FS) from wide (RS) units."""
    if metrics_df is None:
        _no_data(ax, "No metrics CSV"); return
    df = align_metrics(metrics_df, valid)
    snr_m = snr[df["cluster_id"].values]
    sc = ax.scatter(df["duration"], df["halfwidth"],
                    c=snr_m, cmap=CMAP, vmin=SNR_VMIN, vmax=SNR_VMAX,
                    s=8, alpha=0.7, linewidths=0)
    plt.colorbar(sc, ax=ax, label="SNR", pad=0.02)
    ax.set_xlabel("Duration (samples)")
    ax.set_ylabel("Halfwidth (samples)")
    ax.set_title("Duration vs Halfwidth")
    ax.grid(True)
    mask = df["duration"].notna() & df["halfwidth"].notna()
    r = np.corrcoef(df.loc[mask, "duration"], df.loc[mask, "halfwidth"])[0, 1]
    _note(ax, f"r={r:.2f}", "ul")


def panel_pt_ratio(ax, metrics_df, valid):
    """PT ratio histogram — values >1 flag putative inverted / multi-unit waveforms."""
    if metrics_df is None or "pt_ratio" not in metrics_df.columns:
        _no_data(ax, "pt_ratio not found"); return
    df = align_metrics(metrics_df, valid)
    vals = df["pt_ratio"].dropna()
    ax.hist(vals, bins=50, color=PURPLE, alpha=0.85, edgecolor="none")
    med = vals.median()
    ax.axvline(med,  color=GREEN, lw=1.3, ls="--", label=f"median={med:.3f}")
    ax.axvline(1.0,  color=GOLD,  lw=1.0, ls=":",  alpha=0.9, label="PT=1")
    ax.set_xlabel("PT ratio  (peak / |trough|)")
    ax.set_ylabel("# clusters")
    ax.set_title("PT Ratio Distribution")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y")
    inv = (vals > 1).sum()
    _note(ax, f"{inv} units with PT>1  ({100*inv/len(vals):.1f}%)", "ur")


def panel_repol_recovery(ax, fig, metrics_df, snr, valid):
    """Repolarization vs recovery slope, coloured by SNR."""
    if metrics_df is None:
        _no_data(ax, "No metrics CSV"); return
    for col in ("repolarization_slope", "recovery_slope"):
        if col not in metrics_df.columns:
            _no_data(ax, f"Missing column: {col}"); return
    df = align_metrics(metrics_df, valid)
    snr_m = snr[df["cluster_id"].values]
    mask  = df["repolarization_slope"].notna() & df["recovery_slope"].notna()
    sc = ax.scatter(df.loc[mask, "repolarization_slope"],
                    df.loc[mask, "recovery_slope"],
                    c=snr_m[mask], cmap=CMAP, vmin=SNR_VMIN, vmax=SNR_VMAX,
                    s=8, alpha=0.7, linewidths=0)
    plt.colorbar(sc, ax=ax, label="SNR", pad=0.02)
    ax.axhline(0, color=DIM, lw=0.7, ls=":")
    ax.axvline(0, color=DIM, lw=0.7, ls=":")
    ax.set_xlabel("Repolarization slope")
    ax.set_ylabel("Recovery slope")
    ax.set_title("Repolarization vs Recovery Slope")
    ax.grid(True)


def panel_peak_channel(ax, metrics_df, valid, n_channels):
    """
    Histogram of peak_channel indices — shows how units are distributed
    along the probe depth. Gaps can indicate dead/noisy channel banks.
    """
    if metrics_df is None or "peak_channel" not in metrics_df.columns:
        _no_data(ax, "peak_channel not found"); return
    df = align_metrics(metrics_df, valid)
    ax.hist(df["peak_channel"], bins=min(80, n_channels),
            color=GOLD, alpha=0.85, edgecolor="none")
    ax.set_xlabel("Peak channel index")
    ax.set_ylabel("# clusters")
    ax.set_title("Peak Channel Distribution (probe coverage)")
    ax.grid(True, axis="y")
    _note(ax, f"{df['peak_channel'].nunique()} unique channels", "ur")


def panel_trough_peak_offsets(ax, metrics_df, valid):
    """
    Scatter of trough_idx vs peak_idx coloured by duration.
    Well-sorted sessions cluster tightly; scatter indicates temporal jitter
    or heterogeneous waveform shapes.
    """
    if metrics_df is None:
        _no_data(ax, "No metrics CSV"); return
    for col in ("trough_idx", "peak_idx"):
        if col not in metrics_df.columns:
            _no_data(ax, f"Missing column: {col}"); return
    df  = align_metrics(metrics_df, valid)
    dur = df["duration"] if "duration" in df.columns \
          else (df["peak_idx"] - df["trough_idx"]).abs()
    sc  = ax.scatter(df["trough_idx"], df["peak_idx"],
                     c=dur, cmap="plasma", s=8, alpha=0.7, linewidths=0)
    plt.colorbar(sc, ax=ax, label="duration (samples)", pad=0.02)
    lo  = min(df["trough_idx"].min(), df["peak_idx"].min())
    hi  = max(df["trough_idx"].max(), df["peak_idx"].max())
    ax.plot([lo, hi], [lo, hi], color=DIM, lw=0.7, ls="--", alpha=0.6,
            label="peak=trough (ref)")
    ax.set_xlabel("Trough sample index")
    ax.set_ylabel("Peak sample index")
    ax.set_title("Trough vs Peak Sample Position")
    ax.legend(fontsize=7)
    ax.grid(True)


def panel_summary(ax, d, snr, nspk, valid, metrics_df):
    ax.axis("off")
    snr_v  = snr[valid]
    nspk_v = nspk[valid]

    def pct(mask):
        return f"{mask.sum():>4}  ({100*mask.mean():.1f}%)"

    lines = [
        "── Clusters ───────────────────────────",
        f"  Total          : {d['n_clusters']}",
        f"  Valid           : {valid.sum()}",
        f"  Empty (SNR=-1)  : {(~valid).sum()}",
        "",
        "── Waveform geometry ──────────────────",
        f"  Samples/spike   : {d['n_samples']}",
        f"  Channels        : {d['n_channels']}",
        "",
        "── SNR ────────────────────────────────",
        f"  Median          : {np.median(snr_v):.2f}",
        f"  Mean ± SD       : {snr_v.mean():.2f} ± {snr_v.std():.2f}",
        f"  Poor  (<{SNR_THRESHOLDS['poor']})    : {pct(snr_v < SNR_THRESHOLDS['poor'])}",
        f"  Marginal ({SNR_THRESHOLDS['poor']}–{SNR_THRESHOLDS['ok']}) : "
            f"{pct((snr_v >= SNR_THRESHOLDS['poor']) & (snr_v < SNR_THRESHOLDS['ok']))}",
        f"  Good  (≥{SNR_THRESHOLDS['ok']})    : {pct(snr_v >= SNR_THRESHOLDS['ok'])}",
        "",
        "── Spike counts ───────────────────────",
        f"  Median          : {np.median(nspk_v):.0f}",
        f"  Mean ± SD       : {nspk_v.mean():.0f} ± {nspk_v.std():.0f}",
        f"  Min / Max       : {nspk_v.min():.0f} / {nspk_v.max():.0f}",
    ]

    if metrics_df is not None:
        df = align_metrics(metrics_df, valid)
        lines.append("")
        lines.append("── Waveform metrics (medians) ─────────")
        for col, label in [
            ("duration",             "Duration (samp)  "),
            ("halfwidth",            "Halfwidth (samp) "),
            ("pt_ratio",             "PT ratio         "),
            ("repolarization_slope", "Repol. slope     "),
            ("recovery_slope",       "Recovery slope   "),
        ]:
            if col in df.columns:
                med = df[col].dropna().median()
                lines.append(f"  {label}: {med:.4f}")

    ax.text(0.03, 0.97, "\n".join(lines),
            transform=ax.transAxes, va="top", ha="left",
            fontsize=8, fontfamily="monospace", color="#c9d1d9",
            bbox=dict(boxstyle="round,pad=0.55", fc="#161b22",
                      ec="#30363d", alpha=0.95))


# ── main ───────────────────────────────────────────────────────────────────────

def plot_cwave_outputs(data_dir):
    #parser = argparse.ArgumentParser(description="C_Waves diagnostic figure")
    #parser.add_argument("--data_dir",   default=".",
    #                    help="Directory with C_Waves output files (default: .)")
    #parser.add_argument("--out",        default="cwaves_diagnostics.png",
    #                    help="Output path (default: cwaves_diagnostics.png)")
    #parser.add_argument("--max_traces", type=int, default=200,
    #                    help="Max waveforms overlaid in waveform panel (default 200)")
    #args = parser.parse_args()

    data_dir = Path(data_dir)
    d        = load_data(data_dir)
    snr, nspk, valid = build_masks(d)
    metrics  = d["metrics"]

    # ── figure: 4 rows × 3 cols ───────────────────────────────────────────────
    fig = plt.figure(figsize=(21, 22))
    fig.patch.set_facecolor("#0d1117")
    fig.suptitle("C_Waves Diagnostic Summary", fontsize=15,
                 color="#e6edf3", weight="bold", y=0.995)

    gs = gridspec.GridSpec(
        4, 3, figure=fig,
        height_ratios=[1, 1.4, 1, 1],
        hspace=0.45, wspace=0.35,
        left=0.06, right=0.97, top=0.975, bottom=0.04,
    )

    ax_snr_hist  = fig.add_subplot(gs[0, 0])   # R0 — SNR histogram
    ax_spk_hist  = fig.add_subplot(gs[0, 1])   # R0 — spike count histogram
    ax_snr_spk   = fig.add_subplot(gs[0, 2])   # R0 — SNR vs nSpikes

    ax_waveforms = fig.add_subplot(gs[1, :2])  # R1 — waveform overlay (wide)
    ax_summary   = fig.add_subplot(gs[1, 2])   # R1 — text summary

    ax_dur_hw    = fig.add_subplot(gs[2, 0])   # R2 — duration vs halfwidth
    ax_pt        = fig.add_subplot(gs[2, 1])   # R2 — PT ratio histogram
    ax_repol     = fig.add_subplot(gs[2, 2])   # R2 — repol vs recovery

    ax_pk_chan   = fig.add_subplot(gs[3, 0])   # R3 — peak channel distribution
    ax_trpk      = fig.add_subplot(gs[3, 1])   # R3 — trough vs peak sample index
    ax_blank     = fig.add_subplot(gs[3, 2])   # R3 — spare
    ax_blank.axis("off")

    # ── fill panels ───────────────────────────────────────────────────────────
    panel_snr_hist(ax_snr_hist, snr, valid)
    panel_nspikes_hist(ax_spk_hist, nspk, valid)
    panel_snr_vs_nspk(ax_snr_spk, snr, nspk, valid)

    panel_peak_waveforms(ax_waveforms, fig, d["med_wf"], snr, valid,
                         metrics, 200)
    panel_summary(ax_summary, d, snr, nspk, valid, metrics)

    panel_duration_halfwidth(ax_dur_hw, fig, metrics, snr, valid)
    panel_pt_ratio(ax_pt, metrics, valid)
    panel_repol_recovery(ax_repol, fig, metrics, snr, valid)

    panel_peak_channel(ax_pk_chan, metrics, valid, d["n_channels"])
    panel_trough_peak_offsets(ax_trpk, metrics, valid)

    out_path = Path(data_dir / 'cwave_outputs.png')
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"[OK] Figure saved → {out_path.resolve()}")
    plt.show()
    return


if __name__ == "__main__":
    plot_cwave_outputs()
