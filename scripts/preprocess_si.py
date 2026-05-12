"""
preprocess_si.py
Unified ephys preprocessing pipeline: spike sorting + sync processing.

Run order per session:
  1. CatGT          (catgt.do)
  2. OverStrike      (if mouse requires)
  3. DREDge          (motion.do)
  4. Sorter          (sorters.do)
  5. Bombcell        (bombcell.do)
  6. TPrime          (always, requires catgt output)
  7. CWaves          (cwaves.do)
  8. Mean waveform metrics  (mean_waveform_metrics.do)
  9. IBL conversion  (ibl_conversion.do)
"""

import os
import sys
import yaml
import pathlib
import time
import click
from loguru import logger
from pathlib import Path

log_dir = os.environ.get("EPHYS_LOG_DIR", "log")
logger.add(
    f"{log_dir}/preprocess_{{time}}.log",
    colorize=True,
    format="{name} {message}",
    level="INFO",
    rotation="10 MB",
    retention="1 week",
)

from ephys_preprocessing.preprocessing import (
    run_catgt,
    run_overstrike,
    run_dredge,
    run_sorter,
    run_py_bombcell,
    run_tprime,
    run_cwaves,
    run_mean_waveform_metrics,
    run_ibl_ephys_atlas_format,
)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def resolve_processed_dir(input_dir: Path, config: dict) -> Path:
    """
    Derive the processed output directory from the raw input directory,
    mirroring the logic in the original preprocess_spikesort_si.py.
    """
    mouse_name = input_dir.parents[2].name
    if mouse_name == "data":
        mouse_name = input_dir.parents[1].name

    session_name = input_dir.name
    if session_name == "Ephys":
        session_name = input_dir.parents[0].name

    processed_dir = Path(config["output_path"]) / mouse_name / session_name / "Ephys"

    if mouse_name.startswith("MH"):
        processed_dir = Path(str(processed_dir).replace("Axel_Bisi", "Myriam_Hamon"))

    return processed_dir


def find_catgt_folder(ephys_path: Path) -> Path:
    """Return the single catgt_* folder inside an Ephys output directory."""
    catgt_folders = list(ephys_path.glob("catgt_*"))
    if not catgt_folders:
        raise FileNotFoundError(f"No catgt folder found in {ephys_path}")
    if len(catgt_folders) > 1:
        logger.warning(
            f"Multiple catgt folders in {ephys_path}, using: {catgt_folders[0].name}"
        )
    return catgt_folders[0]


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

@logger.catch
def main(input_dir: Path, config_path: Path):
    """
    Run the full pipeline for one session.

    :param input_dir:   Raw session Ephys directory (e.g. /scratch/bisi/data/AB131/Recording/AB131_20240905/Ephys)
    :param config_path: Path to the YAML config file
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    input_dir = Path(input_dir)
    logger.info(f"Starting pipeline for: {input_dir}")
    t0 = time.time()

    def elapsed():
        return time.strftime("%H:%M:%S", time.gmtime(time.time() - t0))

    # Probe count (informational)
    epoch_name = next(d for d in os.listdir(input_dir) if not d.startswith("."))
    n_probes = len([f for f in os.listdir(input_dir / epoch_name) if "imec" in f])
    logger.info(f"Recording with {n_probes} probe(s).")

    # Resolve processed output directory
    processed_dir = resolve_processed_dir(input_dir, config)
    logger.info(f"Output directory: {processed_dir}")
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)

    mouse_name = processed_dir.parents[1].name  # …/{mouse}/{session}/Ephys

    # ------------------------------------------------------------------ #
    # Spike-sorting steps
    # ------------------------------------------------------------------ #

    if config["catgt"]["do"]:
        logger.info("Starting CatGT.")
        run_catgt.main(input_dir, processed_dir, config["catgt"])
        logger.info(f"Finished CatGT in {elapsed()}.")

    # OverStrike — mouse-specific, no config flag needed
    timespans_list = None
    if mouse_name == "PB191":
        timespans_list = [(2350, 2373), (2724, 2778)]
    if timespans_list:
        logger.info("Starting OverStrike.")
        run_overstrike.main(processed_dir, config["overstrike"], timespans_list=timespans_list)
        logger.info(f"Finished OverStrike in {elapsed()}.")

    if config["motion"]["do"]:
        logger.info("Starting DREDge.")
        run_dredge.main(processed_dir, config)
        logger.info(f"Finished DREDge in {elapsed()}.")

    if config["sorters"]["do"]:
        logger.info("Starting sorter.")
        run_sorter.main(processed_dir, config)
        logger.info(f"Finished sorter in {elapsed()}.")

    if config["bombcell"]["do"]:
        logger.info("Starting Bombcell quality metrics.")
        run_py_bombcell.main(processed_dir, config)
        logger.info(f"Finished Bombcell in {elapsed()}.")

    # ------------------------------------------------------------------ #
    # Sync steps — require the catgt folder produced above
    # ------------------------------------------------------------------ #

    catgt_dir = find_catgt_folder(processed_dir)
    logger.info(f"Sync processing from catgt folder: {catgt_dir}")

    logger.info("Starting TPrime.")
    run_tprime.main(catgt_dir, config["tprime"])
    logger.info(f"Finished TPrime in {elapsed()}.")

    if config.get("cwaves", {}).get("do", False):
        logger.info("Starting CWaves.")
        run_cwaves.main(catgt_dir, config["cwaves"])
        logger.info(f"Finished CWaves in {elapsed()}.")

    if config.get("mean_waveform_metrics", {}).get("do", False):
        logger.info("Starting mean waveform metrics.")
        run_mean_waveform_metrics.main(catgt_dir)
        logger.info(f"Finished mean waveform metrics in {elapsed()}.")

    if config.get("ibl_conversion", {}).get("do", False):
        logger.info("Starting IBL conversion.")
        run_ibl_ephys_atlas_format.main(catgt_dir, config)
        logger.info(f"Finished IBL conversion in {elapsed()}.")

    logger.success(f"Pipeline complete in {elapsed()} for {input_dir}.")


# ---------------------------------------------------------------------------
# CLI entry point (used by the SLURM script)
# ---------------------------------------------------------------------------

@click.command()
@click.option(
    "--session",
    required=True,
    help="Absolute or relative session path (e.g. /scratch/bisi/data/AB131/Recording/AB131_.../Ephys)",
)
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to YAML configuration file",
)
@click.option(
    "--data-root",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Prepended to --session when the path is relative",
)
def cli(session, config_path, data_root):
    """
    Run the full ephys preprocessing pipeline (spikesort + sync) for one session.

    Designed to be called once per SLURM array task.
    """
    input_path = Path(session)
    if not input_path.is_absolute():
        if data_root is None:
            with open(config_path) as f:
                cfg = yaml.safe_load(f)
            data_root = Path(cfg["raw_data_path"])
        input_path = data_root / input_path

    main(input_path, config_path)


if __name__ == "__main__":
    cli()
