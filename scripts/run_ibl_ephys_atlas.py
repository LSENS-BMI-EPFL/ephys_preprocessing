#! /usr/bin/env/python3
"""
@author: Jules Lebert
@project: EphysUtils
@file: run_ibl_ephys_atlas.py

Run IBL ephys-atlas formatting on a list of sessions.
CLI version using the same pattern as preprocess_sync_si.py.
"""

import os
import yaml
import time
import click
from datetime import datetime
from joblib import Parallel, delayed
from loguru import logger
from pathlib import Path

# Set up logging — generate filename once so all workers share the same log file
log_dir = os.environ.get('EPHYS_LOG_DIR', 'log')
_log_file = Path(log_dir) / f"run_ibl_ephys_atlas_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
_file_handler_id = logger.add(str(_log_file), colorize=True,
              format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

from ephys_preprocessing.preprocessing import preprocess_ibl_ephys_atlas
from ephys_preprocessing.utils.ephys_utils import transform_input_to_catgt_path, find_catgt_folder


def iter_all_catgt_folders(data_root: Path):
    """
    Iterate over all catgt folders in data_root, structured as:
    data_root/{mouse}/{session}/Ephys/catgt_*

    Skips .DS_Store and other non-directory entries.
    """
    for mouse_dir in sorted(data_root.iterdir()):
        if not mouse_dir.is_dir() or mouse_dir.name == '.DS_Store':
            continue
        for session_dir in sorted(mouse_dir.iterdir()):
            if not session_dir.is_dir() or session_dir.name == '.DS_Store':
                continue
            ephys_dir = session_dir / 'Ephys'
            if not ephys_dir.is_dir():
                continue
            try:
                yield find_catgt_folder(ephys_dir)
            except FileNotFoundError:
                logger.warning(f'No catgt folder found in {ephys_dir}, skipping.')


def main(input_dir, config_file, log_file=None):
    """
    Run IBL ephys-atlas formatting on a single session.
    :param input_dir: path to CatGT preprocessed data (catgt_* folder)
    :param config_file: path to YAML config file
    :param log_file: path to shared log file (used by parallel workers)
    """
    if log_file is not None:
        logger.add(str(log_file), colorize=True, format="{name} {message}",
                   level="INFO", enqueue=True)
    input_dir = Path(input_dir)
    logger.info(f'Processing IBL ephys-atlas format for {input_dir}.')
    start_time = time.time()

    try:
        preprocess_ibl_ephys_atlas.main(input_dir, config_file)
        exec_time = time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time))
        logger.success(f'Finished IBL ephys-atlas formatting in {exec_time} for {input_dir}.')
    except Exception as e:
        logger.error(f'Failed processing {input_dir}: {e}')


@click.command()
@click.option(
    '--input-list',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help='Path to text file containing list of input paths (one per line, same format as spikesort inputs.txt)'
)
@click.option(
    '--config',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help='Path to YAML configuration file'
)
@click.option(
    '--all-sessions',
    is_flag=True,
    help='Process all sessions found in data_root (from config output_path), ignoring --input-list'
)
@click.option(
    '--n-jobs',
    type=int,
    default=1,
    show_default=True,
    help='Number of parallel jobs. Use -1 for all available cores.'
)
def cli(input_list, config, all_sessions, n_jobs):
    """
    IBL ephys-atlas formatting pipeline for electrophysiology data.

    Runs IBL ephys-atlas formatting on a list of sessions. Input paths are
    transformed from spikesort format to find the corresponding catgt folders.

    The YAML config must contain an 'output_path' key pointing to the root
    directory where catgt folders are located (ibl env path).

    Examples:

      \b
      # Process sessions from an input list file
      python run_ibl_ephys_atlas.py --input-list inputs.txt --config config.yaml

      \b
      # Process all sessions found under data_root
      python run_ibl_ephys_atlas.py --all-sessions --config config.yaml

      \b
      # Process all sessions in parallel (4 jobs)
      python run_ibl_ephys_atlas.py --all-sessions --config config.yaml --n-jobs 4
    """
    logger.info('Running IBL ephys-atlas formatting with CLI arguments')

    # Load config to get data root (output_path = catgt parent dir, ibl env path)
    with open(config, 'r') as f:
        config_dict = yaml.load(f, Loader=yaml.FullLoader)

    data_root = Path(config_dict['output_path'])
    logger.info(f'Looking for catgt folders in: {data_root}')

    if all_sessions:
        logger.info('Running in all-sessions mode: iterating over all mice/sessions in data_root')
        catgt_paths = list(iter_all_catgt_folders(data_root))
    else:
        if not input_list:
            raise click.UsageError('--input-list is required unless using --all-sessions')

        with open(input_list, 'r') as f:
            input_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        logger.info(f'Loaded {len(input_paths)} input paths from {input_list}')

        catgt_paths = []
        for input_path in input_paths:
            try:
                catgt_path = transform_input_to_catgt_path(input_path, data_root)
                logger.info(f'Transformed path: {input_path} -> {catgt_path}')
                catgt_paths.append(catgt_path)
            except FileNotFoundError as e:
                logger.error(f'Skipping {input_path}: {e}')

    logger.info(f'Processing {len(catgt_paths)} session(s) with n_jobs={n_jobs}')
    # Remove all loguru sinks before spawning workers — file handles can't be pickled.
    # This includes sinks added by imported modules (e.g. preprocess_ibl_ephys_atlas).
    # Worker processes re-add their own sinks on module import.
    logger.remove()
    # Parallel(n_jobs=n_jobs)(delayed(main)(p, config, _log_file) for p in catgt_paths)
    for p in catgt_paths:
        main(p, config, _log_file)


if __name__ == '__main__':
    cli()
