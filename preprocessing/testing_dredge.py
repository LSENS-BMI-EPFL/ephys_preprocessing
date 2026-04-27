"""
DREDge Motion Estimation Pipeline
==================================
Takes a preprocessed .ap.bin (or .lf.bin) SpikeGLX file, estimates motion
with DREDge, and saves the Motion object + a drift figure.

Nothing is written to the traces — no corrected binary, no sorting.
The saved Motion object can later be loaded and passed to
si.interpolate_motion() before spike sorting.
"""

import gc
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import psutil

# ---------------------------------------------------------------------------
# Root logger — console only at startup; per-job file handlers added later
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# JOB_KWARGS
# n_jobs=4 instead of 8:  each worker loads ~1–2 GB of chunk data;
# 8 workers on a 45-min recording was likely hitting 16–20 GB and triggering
# the Windows OOM killer, which causes BrokenProcessPool.
# Lower this further (n_jobs=2) if you still see crashes.
# ---------------------------------------------------------------------------
JOB_KWARGS = dict(
    n_jobs=24,
    chunk_duration="1s",
    progress_bar=True,
)


# ---------------------------------------------------------------------------
# Per-job log file helper
# ---------------------------------------------------------------------------
def _make_job_log_handler(log_path: Path) -> logging.FileHandler:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S",
        )
    )
    return fh


def _close_handler_safe(fh: logging.FileHandler):
    """Close a FileHandler without raising even if the underlying stream
    was corrupted by a dying subprocess pool (OSError EINVAL on Windows)."""
    root_logger = logging.getLogger()
    try:
        root_logger.removeHandler(fh)
    except Exception:
        pass
    try:
        fh.flush()
    except Exception:
        pass
    try:
        fh.stream.close()
    except Exception:
        pass
    try:
        fh.close()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# SpikeGLX binary reader
# ---------------------------------------------------------------------------
def _read_spikeglx_bin(bin_file: Path, load_sync_channel: bool = False):
    import spikeinterface.full as si
    from probeinterface import read_spikeglx as pi_read_spikeglx

    meta_file = bin_file.with_suffix(".meta")
    if not meta_file.exists():
        raise FileNotFoundError(f"Meta file not found: {meta_file}")

    meta = {}
    with open(meta_file, "r") as f:
        for line in f:
            line = line.strip()
            if "=" in line:
                k, v = line.split("=", 1)
                meta[k.strip()] = v.strip()

    sampling_frequency = float(meta["imSampRate"])
    n_channels_total   = int(meta["nSavedChans"])
    n_channels         = n_channels_total if load_sync_channel else n_channels_total - 1

    rec = si.read_binary(
        bin_file,
        sampling_frequency=sampling_frequency,
        num_channels=n_channels,
        dtype="int16",
        gain_to_uV=1.0,
        offset_to_uV=0.0,
        channel_ids=[str(i) for i in range(n_channels)],
    )

    try:
        probe = pi_read_spikeglx(meta_file)
        rec   = rec.set_probe(probe)
        log.info("Probe attached from meta file.")
    except Exception as e:
        log.warning("Could not attach probe: %s — depth axis may be missing.", e)

    return rec


# ---------------------------------------------------------------------------
# Memory helpers
# ---------------------------------------------------------------------------
def _log_ram(label: str = ""):
    proc     = psutil.Process(os.getpid())
    mb       = proc.memory_info().rss / 1024**2
    total_gb = psutil.virtual_memory().total / 1024**3
    used_pct = psutil.virtual_memory().percent
    log.info(
        "RAM %s: %.0f MB process | %.1f%% of %.0f GB system used",
        label, mb, used_pct, total_gb,
    )


def _log_vram():
    try:
        import torch
        if torch.cuda.is_available():
            alloc    = torch.cuda.memory_allocated() / 1024**2
            reserved = torch.cuda.memory_reserved()  / 1024**2
            log.info("VRAM: %.0f MB allocated / %.0f MB reserved", alloc, reserved)
    except ImportError:
        pass


def _clear_gpu():
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
            log.info("GPU cache cleared.")
            _log_vram()
    except ImportError:
        pass


def _wait_for_ram(threshold_gb: float = 12.0, poll_interval: float = 30.0):
    """Block until at least threshold_gb of system RAM is free."""
    while True:
        vm      = psutil.virtual_memory()
        free_gb = vm.available / 1024**3
        if free_gb >= threshold_gb:
            return
        log.warning(
            "Only %.1f / %.1f GB free — waiting %.0fs …",
            free_gb, vm.total / 1024**3, poll_interval,
        )
        time.sleep(poll_interval)


# ---------------------------------------------------------------------------
# Core run() — one probe
# ---------------------------------------------------------------------------
def run(
    bin_file: Path,
    output_folder: Path,
    preset: str = "dredge",
    use_lfp: bool = False,
    overwrite: bool = True,
):
    import spikeinterface.full as si
    from spikeinterface.sortingcomponents.motion import estimate_motion

    bin_file      = Path(bin_file)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    motion_folder  = output_folder / "motion"
    figures_folder = output_folder / "figures"
    figures_folder.mkdir(exist_ok=True)

    if overwrite and motion_folder.exists():
        shutil.rmtree(motion_folder)

    log.info("Loading %s", bin_file.name)
    rec = _read_spikeglx_bin(bin_file, load_sync_channel=False) #stream name?
    log.info("Recording: %s", rec)
    _log_ram("after loading recording")

    if use_lfp:
        log.info("Estimating motion with dredge_lfp …")
        rec_lfp = si.bandpass_filter(rec, freq_min=0.5, freq_max=500.0)
        rec_lfp = si.phase_shift(rec_lfp)
        rec_lfp = si.resample(rec_lfp, resample_rate=250)
        del rec

        motion = estimate_motion(
            rec_lfp,
            method="dredge_lfp",
            rigid=True,
            progress_bar=True,
        )
        del rec_lfp

    else:
        log.info(
            "Estimating motion with dredge_ap (preset='%s', n_jobs=%d) …",
            preset, JOB_KWARGS["n_jobs"],
        )

        NT = 64 * 1024 + 64
        FS = rec.sampling_frequency
        rec.reset_times()
        rec = si.highpass_filter(recording=rec, freq_min=300.)

        motion_dict = {
            'bin_s': NT / FS,
            'histogram_time_smooth_s': NT / FS
        }

        ## Mayo
        #rec_corrected, motion = si.correct_motion(recording=rec,
        #                                          preset=preset,
        #                                          folder=motion_folder,
        #                                          output_motion_info=True,
        #                                          n_jobs=24,
        #                                          progress_bar=True,
        #                                          estimate_motion_kwargs=motion_dict,
        #                                          overwrite=overwrite)
        #print(type(rec_corrected), type(motion))
        #_save_motion_figure(motion, figures_folder)

        ## Myriam
        motion, motion_info = si.compute_motion(
            rec,
            preset=preset,
            folder=motion_folder,
            output_motion_info=True,
            overwrite=overwrite,
            estimate_motion_kwargs=motion_dict,
            **JOB_KWARGS,
        )

        del rec
        del motion_info

    log.info("Motion estimated: %s", motion)
    _log_ram("after motion estimation")
    _log_vram()

    _save_motion_figure(motion, figures_folder)
    del motion

    gc.collect()
    _clear_gpu()
    _log_ram("after full cleanup")
    log.info("Done. Motion saved in %s", motion_folder)


# ---------------------------------------------------------------------------
# Motion figure
# ---------------------------------------------------------------------------
def _save_motion_figure(motion, figures_folder: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 4))
        disp  = motion.displacement[0]
        times = motion.temporal_bins_s[0]

        if disp.ndim == 2 and disp.shape[1] > 1:
            for i, d in enumerate(motion.spatial_bins_um):
                ax.plot(times, disp[:, i], lw=0.8, alpha=0.7, label=f"{d:.0f} µm")
            ax.legend(fontsize=6, ncol=4, loc="upper right")
        else:
            ax.plot(times, disp if disp.ndim == 1 else disp[:, 0], lw=1)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Displacement (µm)")
        ax.set_title("DREDge estimated motion")
        fig.tight_layout()

        out = figures_folder / "dredge_motion.png"
        fig.savefig(out, dpi=150)
        log.info("Figure saved → %s", out)
        plt.close(fig)
        plt.close("all")

    except Exception as e:
        log.warning("Could not save motion figure: %s", e)


# ---------------------------------------------------------------------------
# Directory crawl
# ---------------------------------------------------------------------------
def collect_all_jobs(
    data_root: Path,
    use_lfp: bool = False,
):
    data_root = Path(data_root)
    pattern   = "*corrected*.lf.bin" if use_lfp else "*corrected*.ap.bin"
    jobs      = []

    for mouse_dir in sorted(data_root.iterdir())[20:30]:
        if not mouse_dir.is_dir():
            continue
        for session_dir in sorted(mouse_dir.iterdir()):
            ephys_dir = session_dir / "Ephys"
            if not ephys_dir.is_dir():
                continue
            for catgt_dir in sorted(ephys_dir.glob("catgt_*")):
                if not catgt_dir.is_dir():
                    continue
                for imec_dir in sorted(catgt_dir.glob("*imec*")):
                    if not imec_dir.is_dir():
                        continue
                    bin_files = [f for f in imec_dir.glob(pattern) if f.is_file()]
                    if not bin_files:
                        continue
                    if len(bin_files) > 1:
                        log.warning(
                            "Multiple matches in %s — using: %s",
                            imec_dir.name, bin_files[0].name,
                        )
                    jobs.append((bin_files[0], imec_dir / "dredge"))
                    log.info("Queued: %s", bin_files[0].relative_to(data_root))

    return jobs


# ---------------------------------------------------------------------------
# run_all
# ---------------------------------------------------------------------------
def run_all(
    data_root: Path,
    preset: str = "dredge",
    use_lfp: bool = False,
    overwrite: bool = False,
    skip_existing: bool = True,
    ram_threshold_gb: float = 12.0,
):
    data_root = Path(data_root)
    jobs = collect_all_jobs(data_root, use_lfp=use_lfp)

    if not jobs:
        log.warning("No corrected .bin files found under %s", data_root)
        return

    if skip_existing and not overwrite:
        before = len(jobs)
        jobs   = [(b, o) for b, o in jobs if not (o / "motion").exists()]
        if skipped := before - len(jobs):
            log.info("Skipping %d already-processed probe(s).", skipped)

    log.info(
        "Running DREDge on %d probe(s) sequentially (n_jobs=%d per probe).",
        len(jobs), JOB_KWARGS["n_jobs"],
    )

    summary: list[dict] = []

    for job_idx, (bin_file, output_folder) in enumerate(jobs, start=1):
        label = bin_file.parent.name

        # Open per-job log AFTER previous job's subprocess pool is fully gone
        log_path = output_folder / "dredge_motion.log"
        job_fh   = _make_job_log_handler(log_path)
        logging.getLogger().addHandler(job_fh)

        t0    = time.time()
        error = None

        try:
            log.info("=" * 60)
            log.info("JOB %d / %d  —  %s", job_idx, len(jobs), label)
            log.info("bin  : %s", bin_file)
            log.info("out  : %s", output_folder)
            log.info("n_jobs (workers): %d", JOB_KWARGS["n_jobs"])
            log.info("=" * 60)

            _wait_for_ram(threshold_gb=ram_threshold_gb)
            _log_ram(f"[{label}] before run")

            run(
                bin_file=bin_file,
                output_folder=output_folder,
                preset=preset,
                use_lfp=use_lfp,
                overwrite=overwrite,
            )

            elapsed = time.time() - t0
            log.info("[%s] Completed in %.1f min.", label, elapsed / 60)

        except Exception as exc:
            elapsed = time.time() - t0
            error   = exc
            try:
                log.error(
                    "[%s] FAILED after %.1f min: %s",
                    label, elapsed / 60, exc,
                    exc_info=True,
                )
            except Exception:
                pass  # handler stream may be broken — don't crash the loop

        finally:
            gc.collect()
            _clear_gpu()
            try:
                _log_ram(f"[{label}] after cleanup")
            except Exception:
                pass
            _close_handler_safe(job_fh)  # safe even if subprocess pool died

        summary.append(
            {
                "label":       label,
                "status":      "FAILED" if error else "OK",
                "elapsed_min": round(elapsed / 60, 1),
                "log":         str(log_path),
                "error":       str(error) if error else "",
            }
        )

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    log.info("")
    log.info("=" * 60)
    log.info("SUMMARY  (%d jobs)", len(summary))
    log.info("=" * 60)
    for s in summary:
        suffix = f"  ← {s['error']}" if s["error"] else ""
        log.info("  [%s] %s  (%.1f min)%s", s["status"], s["label"], s["elapsed_min"], suffix)
    log.info("Logs written alongside each probe's dredge/ folder.")
    log.info("All done.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    data_root = r"M:\analysis\Axel_Bisi\data"
    #
    # run_all(
    #     data_root=data_root,
    #     preset="dredge",
    #     use_lfp=False,
    #     overwrite=False,
    #     skip_existing=True,
    #     ram_threshold_gb=12.0,
    # )
    bin_file = r"M:\analysis\Axel_Bisi\data\AB133\AB133_20241105_111234\Ephys\catgt_AB133_g0\AB133_g0_imec2\AB133_g0_tcat_corrected.imec2.ap.bin"
    output_folder = r"M:\analysis\Myriam_Hamon\Presentations\IBL_meetings\test\dredge\AB133_imec2"
    run(
        bin_file=bin_file,
        output_folder=output_folder,
        preset="dredge",
        use_lfp=False,
        overwrite=False,
    )