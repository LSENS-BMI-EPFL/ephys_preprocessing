import os
import shutil
from pathlib import Path
from loguru import logger
import numpy as np
import yaml
import spikeinterface.extractors as se
import spikeinterface.full as si

import spikeinterface.preprocessing as spre
from spikeinterface.exporters import export_to_ibl_gui, export_to_phy
from spikeinterface.core import ChannelSparsity
from spikeinterface.extractors import read_kilosort_as_analyzer

from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording, get_exp_datetime
from ephys_preprocessing.preprocessing.get_artifact_times import run_tprime_alignment
from ephys_preprocessing.preprocessing.spikeinterface_modules.artifact_correction import artifact_correction

def main(input_dir, config):
    # epoch_name = os.listdir(input_dir)[0]
    input_dir = Path(input_dir)
    date = get_exp_datetime(input_dir)
    epoch_name = list(input_dir.glob('catgt*'))[0].name
    probe_folders = (input_dir / epoch_name).glob("*imec*")
    probe_folders = [str(p) for p in probe_folders]
    logger.info('Data to spike-sort: {}'.format(probe_folders))
    n_probes = len(probe_folders)

    sorter_configs = config['sorters']

    if 'kilosort2' in sorter_configs['sorter_list']:
        si.Kilosort2Sorter.set_kilosort2_path(sorter_configs['kilosort2']['sorter_path'])

    for probe_id in range(n_probes):
        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[1]
        if not check_if_valid_recording(config, mouse_id, probe_id, date):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)
        preprocessed_path = Path(probe_path) / 'preprocess'
        overwrite_preprocessed = config['sorters']['overwrite_preprocessed']

        if preprocessed_path.exists() and any(preprocessed_path.iterdir()) and not overwrite_preprocessed:
            recording = si.load_extractor(preprocessed_path) #TODO: this throws an error
        else:
            # Get recording
            recording = se.read_spikeglx(probe_path, stream_id=f'imec{probe_id}.ap')

            # Preprocessing steps before spike sorting
            # --------------------------------------
            if config['artifact_correction']['do']:
                logger.info('Correcting magnetic artifact before sorting...')

                artifact_times = run_tprime_alignment(
                    input_dir=input_dir,
                    probe_id=probe_id,
                    force_rerun=True,
                    tprime_config=config['tprime'],
                )
                #recording = artifact_correction(
                #    recording=recording,
                #    artifact_times=artifact_times,
                #    window_ms=config['artifact_correction']['window_ms'],
                #)

                artifact_times_frame = artifact_times * recording.sampling_frequency
                recording = si.remove_artifacts(
                    recording, 
                    list_triggers=artifact_times_frame,
                    ms_before=config['artifact_correction']['ms_before'],
                    ms_after=config['artifact_correction']['ms_after'],
                    mode=config['artifact_correction']['mode'],
                )

            # Drift correction before sorting to avoid problem on the clusters
            # where drift correction hangs forever (kilosort4)
            if config['drift_correction']['do']:
                logger.info('Correcting drift before sorting...')
                recording = si.correct_motion(
                    recording=recording,
                    preset=config['drift_correction']['preset'],
                    **config['sorters']['job_kwargs']
                )

            logger.info('Saving artifact/drift preprocessed recording raw data')
            recording = recording.save(
                folder = preprocessed_path,
                format='binary',
                overwrite=True,
                **config['sorters']['job_kwargs'],
            )


        # Run spike sorting
        for sorter in sorter_configs['sorter_list']:


            folder = Path(probe_path) / sorter

            logger.info(f'Running {sorter} for IMEC probe {probe_id}')
            sorting = si.run_sorter(
                sorter_name=sorter,
                recording=recording,
                remove_existing_folder=True,
                folder=folder,
                verbose=True,
                singularity_image=sorter_configs[sorter]["singularity_image"],
                **sorter_configs[sorter]['sorter_params'],
            )

if __name__ == "__main__":
    input_dir = Path('/Volumes/Petersen-Lab/analysis/Axel_Bisi/data/PB191/PB191_20241210_110601/Ephys/')
    config_file = '/Users/lebert/home/code/preprocessing_tools/ephys_preprocessing/ephys_preprocessing/preprocessing/preprocess_config_si.yaml'
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)


    raw_path = Path('/Volumes/Petersen-Lab/data/PB191/Recording/Ephys/PB191_20241210_110601/PB191_g0')
    stream_names = stream_names, stream_ids = si.get_neo_streams('spikeglx', raw_path)
    rec = si.read_spikeglx(raw_path, "imec2.ap")
    print(rec)

    main(input_dir, config)