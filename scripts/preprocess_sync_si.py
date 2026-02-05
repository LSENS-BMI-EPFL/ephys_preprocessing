#! /usr/bin/env/python3
"""
@author: Axel Bisi, Jules Lebert
@project: EphysUtils
@file: preprocess_sync_si.py
@time: 8/2/2023 6:28 PM

Run TPrime, Cwaves, and waveform metrics on preprocessed spike data.
CLI version for Docker/Singularity deployment.
"""

# Imports
import os
import yaml
import time
import click
from loguru import logger
from pathlib import Path

# Set up logging - use environment variable for log directory if available
log_dir = os.environ.get('EPHYS_LOG_DIR', 'log')
logger.add(f"{log_dir}/preprocess_sync_{{time}}.log", colorize=True,
              format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
from ephys_preprocessing.preprocessing import (
     run_tprime,
     run_cwaves,
     run_mean_waveform_metrics,
    #  run_lfp_analysis,
)


def find_catgt_folder(ephys_path: Path) -> Path:
    """
    Find the catgt folder inside an Ephys directory.

    The catgt folder follows the pattern: catgt_{mouse_id}_{date}*_g*

    :param ephys_path: Path to the Ephys directory (processed output)
    :return: Path to the catgt folder
    :raises FileNotFoundError: If no catgt folder is found
    """
    catgt_folders = list(ephys_path.glob('catgt_*'))

    if not catgt_folders:
        raise FileNotFoundError(f"No catgt folder found in {ephys_path}")

    if len(catgt_folders) > 1:
        logger.warning(f"Multiple catgt folders found in {ephys_path}, using first one: {catgt_folders[0].name}")

    return catgt_folders[0]


def transform_input_to_catgt_path(input_path: str, data_root: Path) -> Path:
    """
    Transform an input path from inputs.txt format to the catgt folder path.

    Handles two input formats:
    1. Normal format: JL005/Recording/JL005_20250520_142542/Ephys
       -> JL005/JL005_20250520_142542/Ephys/catgt_*

    2. PB format: PB197/Recording/Ephys/PB197_20241216_155436
       -> PB197/PB197_20241216_155436/Ephys/catgt_*

    :param input_path: Path from inputs.txt
    :param data_root: Root directory for processed data (e.g., /scratch/lebert/ephys_output)
    :return: Full path to the catgt folder
    """
    parts = input_path.strip('/').split('/')

    # Detect format based on path structure
    # PB format: {mouse}/Recording/Ephys/{session} (Ephys is 3rd element, session is 4th)
    # Normal format: {mouse}/Recording/{session}/Ephys (session is 3rd element, Ephys is 4th)

    if len(parts) >= 4 and parts[2] == 'Ephys':
        # PB format: PB197/Recording/Ephys/PB197_20241216_155436
        mouse_name = parts[0]
        session_name = parts[3]
        processed_path = f"{mouse_name}/{session_name}/Ephys"
    else:
        # Normal format: JL005/Recording/JL005_20250520_142542/Ephys
        # Just remove '/Recording/' from the path
        processed_path = input_path.replace('/Recording/', '/')

    # Build full path to Ephys directory
    ephys_path = data_root / processed_path

    # Find the catgt folder
    catgt_path = find_catgt_folder(ephys_path)

    return catgt_path


@logger.catch
def main(input_dir, config):
    """
    Run TPrime and Cwaves on preprocessed spike data.
    :param input_dir: path to CatGT preprocessed data (catgt_* folder)
    :param config: config dictionary or path to config file
    :return:
    """
    input_dir = Path(input_dir)

    # Load config if it's a path
    if isinstance(config, (str, Path)):
        with open(config, 'r') as f:
            config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Processing sync data from {}.'.format(input_dir))
    start_time = time.time()

    # Run TPrime
    logger.info('Starting Tprime.')
    run_tprime.main(input_dir, config['tprime'])
    logger.info('Finished Tprime in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Run Cwaves
    logger.info('Starting Cwaves.')
    run_cwaves.main(input_dir, config['cwaves'])
    logger.info('Finished Cwaves in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Run mean waveform metrics
    logger.info('Starting mean waveform metrics.')
    run_mean_waveform_metrics.main(input_dir)
    logger.info('Finished mean waveform metrics in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # LFP analysis for depth estimation (optional, currently disabled)
    # logger.info('Starting LFP analysis.')
    # run_lfp_analysis.main(input_dir)
    # logger.info('Finished LFP analysis in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    exec_time_hhmmss = time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))
    logger.success(f'Finished sync processing in {exec_time_hhmmss} for {input_dir}.')

    return


@click.command()
@click.option(
    '--input-list',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to text file containing list of input paths (one per line, same format as spikesort inputs.txt)'
)
@click.option(
    '--config',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to YAML configuration file'
)
@click.option(
    '--legacy-mode',
    is_flag=True,
    help='Use hardcoded paths from script (for local debugging)'
)
def cli(input_list, config, legacy_mode):
    """
    Sync processing pipeline for electrophysiology data.

    Runs TPrime, Cwaves, and waveform metrics on CatGT-processed data.

    This script transforms input paths from the spikesort format to find
    the corresponding catgt folders in the processed output directory.

    Examples:

      \b
      # Process sessions from an input list file
      python preprocess_sync_si.py --input-list inputs.txt --config config.yaml

      \b
      # Use legacy mode with hardcoded paths (for debugging)
      python preprocess_sync_si.py --legacy-mode

      \b
      # Docker/Singularity usage
      singularity exec ephys-pipeline.sif python preprocess_sync_si.py \\
        --input-list /mnt/config/inputs.txt --config /mnt/config/config.yaml
    """

    # Legacy mode with hardcoded paths
    if legacy_mode:
        logger.info('Running in legacy mode with hardcoded paths')
        data_path = Path('/home/lebert/lsens_srv/analysis/Jules_Lebert/data')
        input_paths = [
            # Already in catgt format for legacy mode
            'JL005/JL005_20250520_142542/Ephys/catgt_JL005_20250520_g0',
            'JL002/JL002_20250522_111333/Ephys/catgt_JL002_20250522_2_g0',
            'JL002/JL002_20250523_101907/Ephys/catgt_JL002_20250523_g0',
            'JL007/JL007_20250523_144849/Ephys/catgt_JL007_20250523_g0',
            'JL006/JL006_20250601_104051/Ephys/catgt_JL006_20250601_g0',
            'JL006/JL006_20250602_122916/Ephys/catgt_JL006_20250602_g0',
            'JL007/JL007_20250603_150143/Ephys/catgt_JL007_20250603_2_g0',
            'JL007/JL007_20250605_145217/Ephys/catgt_JL007_20250605_2_g0'
        ]
        config_path = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/scripts/preprocess_config_si.yaml')

        for input_path in input_paths:
            main(data_path / input_path, config_path)
    else:
        # New mode with CLI arguments
        if not input_list or not config:
            raise click.UsageError('--input-list and --config are required unless using --legacy-mode')

        logger.info('Running with CLI arguments')

        # Load config to get data paths
        with open(config, 'r') as f:
            config_dict = yaml.load(f, Loader=yaml.FullLoader)

        # Use output_path as the data root (where catgt folders are)
        data_root = Path(config_dict['output_path'])

        # Read input list from file
        with open(input_list, 'r') as f:
            input_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        logger.info(f'Loaded {len(input_paths)} input paths from {input_list}')
        logger.info(f'Looking for catgt folders in: {data_root}')

        # Process each input
        for input_path in input_paths:
            try:
                # Transform the input path to find the catgt folder
                catgt_path = transform_input_to_catgt_path(input_path, data_root)
                logger.info(f'Transformed path: {input_path} -> {catgt_path}')
                main(catgt_path, config_dict)
            except FileNotFoundError as e:
                logger.error(f'Skipping {input_path}: {e}')
                continue


if __name__ == '__main__':
    cli()
