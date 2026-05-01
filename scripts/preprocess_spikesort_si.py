# Imports
import os
import yaml
import pathlib
import time
import click
from loguru import logger
from pathlib import Path
logger.add("log/preprocess_spikesort_{time}.log", colorize=True,
           format="{name} {message}", level="INFO", rotation="10 MB", retention="1 week")

# Import submodules
from ephys_preprocessing.preprocessing import (
    run_catgt,
    run_overstrike,
    run_dredge,
    run_sorter,
    run_py_bombcell,
)

@logger.catch
def main(input_dir, config_file):
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    logger.info('Preprocessing data from {}.'.format(input_dir))
    start_time = time.time()

    # Get epoch number and run name
    epoch_name = [dir for dir in os.listdir(input_dir) if not dir.startswith('.')][0]
    # Find only directory names with "imec" in it
    n_probes = len([f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f])
    logger.info('Recording using {} probe(s).'.format(n_probes))

    # Create output folder
    #mouse_name = input_dir.split('\\')[2]
    #session_name = input_dir.split('\\')[-2]
    mouse_name = input_dir.parents[2].name
    if mouse_name == 'data':
        mouse_name = input_dir.parents[1].name

    session_name = input_dir.name
    if session_name == 'Ephys':
        session_name = input_dir.parents[0].name
    print('Mouse name', mouse_name)
    print('Session name', session_name)
    processed_dir = os.path.join(config['output_path'], mouse_name, session_name, 'Ephys')
    if mouse_name.startswith('MH'):
        processed_dir = processed_dir.replace('Axel_Bisi', 'Myriam_Hamon')

    logger.info('Processed data will be saved to {}.'.format(processed_dir))
    pathlib.Path(processed_dir).mkdir(parents=True, exist_ok=True)


    # Run CatGT
    if config['catgt']['do']:
        logger.info('Starting CatGT.')
        run_catgt.main(input_dir, processed_dir, config['catgt'])
        logger.info('Finished CatGT in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # Optionally, run OverStrike
    timespans_list = None
    if mouse_name == 'PB191':
        timespans_list = [(2350, 2373), (2724, 2778)] #add timespans to overtrike config here (trials will be ignored at NWB conversion)

    if timespans_list:
        logger.info('Starting OverStrike.')
        run_overstrike.main(processed_dir, config['overstrike'], timespans_list=timespans_list)
        logger.info('Finished OverStrike in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))

    # # Run DREDge
    if config['motion']['do']:
        logger.info('Starting DREDge.')
        run_dredge.main(processed_dir, config)
        logger.info("Finished DREDge in {}.".format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))


    # # Run Kilosort
    if config['sorters']['do']:
        logger.info('Starting Kilosort.')
        run_sorter.main(processed_dir, config)
        logger.info("Finished Kilosort in {}.".format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))


    # # Run quality metrics e.g. bombcell
    if config['bombcell']['do']:
        logger.info('Starting bombcell quality metrics.')
        run_py_bombcell.main(processed_dir, config)
        logger.info('Finished bombcell quality metrics in {}.'.format(time.strftime('%H:%M:%S', time.gmtime(time.time()-start_time))))


@click.command()
@click.option(
    '--input-list',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to text file containing list of input paths (one per line, relative to data root)'
)
@click.option(
    '--config',
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help='Path to YAML configuration file'
)
# @click.option(
#     '--data-root',
#     type=click.Path(exists=True, file_okay=False, path_type=Path),
#     default='/mnt/lsens/data',
#     help='Root directory for input data (default: /mnt/lsens/data)',
#     show_default=True
# )
@click.option(
    '--legacy-mode',
    is_flag=True,
    help='Use hardcoded paths from script (for local debugging)'
)
def cli(input_list, config, legacy_mode):
    """
    Electrophysiology preprocessing pipeline with SpikeInterface.

    Runs CatGT, OverStrike (optional), Kilosort, and Bombcell quality metrics.

    Examples:

      \b
      # Process sessions from an input list file
      python preprocess_spikesort_si.py --input-list inputs.txt --config config.yaml

      \b
      # Use legacy mode with hardcoded paths (for debugging)
      python preprocess_spikesort_si.py --legacy-mode

      \b
      # Docker usage
      docker run -v /data:/mnt/data -v /output:/mnt/output ephys-pipeline \\
        --input-list /mnt/data/inputs.txt --config /mnt/data/config.yaml
    """

    # Legacy mode with hardcoded paths
    if legacy_mode:
        logger.info('Running in legacy mode with hardcoded paths')

        from pathlib import Path
        from concurrent.futures import ThreadPoolExecutor

        analysis_path = Path(r"M:\analysis\Axel_Bisi\data")
        data_path = Path(r"M:\data")

        def process_session(session_dir):
            ephys = session_dir / "Ephys"
            if not ephys.exists():
                return None
            if any(ephys.glob("catgt_*/*/*tcat_corrected*")):
                mouse = session_dir.parent.name
                session = session_dir.name
                return data_path / mouse / "Recording" / session / "Ephys"

        sessions = [
            session
            for mouse in analysis_path.iterdir() if mouse.is_dir() and 'AB' in mouse.name
            for session in mouse.iterdir() if session.is_dir()
        ]

        with ThreadPoolExecutor() as ex:
            input_paths = sorted(p for p in ex.map(process_session, sessions) if p)

        for p in input_paths:
            print(p)

       #input_paths = []
       #for f in sorted(path_analysis.rglob("*tcat_corrected*")):
       #    ephys = next((p for p in f.parents if p.name == "Ephys"), None)
       #    if ephys is None:
       #        continue
       #    new_path = data_path / ephys.relative_to(path_analysis)
       #    if new_path not in input_paths and ('AB' in str(new_path)):  # or 'MH' in str(new_path)
       #        print(new_path)
       #        input_paths.append(new_path)

        #input_paths = [
            #'MH035/Recording/MH035_20250514_142713/Ephys',
            #'MH031/Recording/MH031_20250507_104425/Ephys',
            #'AB131/Recording/AB131_20240905_123601/Ephys',
            # 'PB191/Recording/Ephys/PB191_20241210_110601',
            # 'PB192/Recording/Ephys/PB192_20241211_113347',
            # 'PB201/Recording/Ephys/PB201_20241212_192123',
            # 'PB193/Recording/Ephys/PB193_20241218_135125',
            # 'PB194/Recording/Ephys/PB194_20241218_161235',
            # 'PB195/Recording/Ephys/PB195_20241214_114010',
            # 'PB196/Recording/Ephys/PB196_20241217_144715',
            # 'PB197/Recording/Ephys/PB197_20241216_155436',
            # 'PB198/Recording/Ephys/PB198_20241213_142448',
            # 'PB200/Recording/Ephys/PB200_20241216_112934',
            # 'RD076/Recording/RD076_20250214_125235/Ephys'
            # 'RD077/Recording/RD077_20250219_183425/Ephys',
            # 'RD077/Recording/RD077_20250221_102024/Ephys',
            # 'RD072/Recording/RD072_20250305_131521/Ephys',
            # 'JL002/Recording/JL002_20250507_135553/Ephys',
            # 'JL005/Recording/JL005_20250520_142542/Ephys/',
            # 'JL002/Recording/JL002_20250522_111333/Ephys/',
            # 'JL002/Recording/JL002_20250523_101907/Ephys/',
            # 'JL007/Recording/JL007_20250523_144849/Ephys/',
            # 'JL006/Recording/JL006_20250601_104051/Ephys/',
            # 'JL006/Recording/JL006_20250602_122916/Ephys/',
            #'JL007/Recording/JL007_20250603_150143/Ephys/',
            #'JL007/Recording/JL007_20250605_145217/Ephys/'
        #]
        print(input_paths)
        #config_path = Path('/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/scripts/preprocess_config_si.yaml')
        config_path = Path(r"C:\Users\bisi\Github\ephys_preprocessing\scripts\preprocess_config_si.yaml")

    else:
        # New mode with CLI arguments
        if not input_list or not config:
            raise click.UsageError('--input-list and --config are required unless using --legacy-mode')

        logger.info('Running with CLI arguments')
        config_path = config
        with open(config_path, 'r') as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        data_path = Path(config['raw_data_path'])

        # Read input list from file
        with open(input_list, 'r') as f:
            input_paths = [line.strip() for line in f if line.strip() and not line.startswith('#')]

        logger.info(f'Loaded {len(input_paths)} input paths from {input_list}')

    # Process each input
    for input_path in input_paths:
        print(data_path / input_path, config_path)
        main(data_path / input_path, config_path)


if __name__ == '__main__':
    cli()