import os
import shutil
from pathlib import Path
from loguru import logger
import yaml
import spikeinterface.extractors as se
import spikeinterface.full as si

from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording
from ephys_preprocessing.preprocessing.get_artifact_times import run_tprime_alignment
from ephys_preprocessing.preprocessing.spikeinterface_modules.artifact_correction import artifact_correction

def main(input_dir, config):
    # epoch_name = os.listdir(input_dir)[0]
    input_dir = Path(input_dir)
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
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)
        preprocessed_path = Path(probe_path) / 'preprocess'
        overwrite_preprocessed = config['sorters']['overwrite_preprocessed']

        if preprocessed_path.exists() and any(preprocessed_path.iterdir()) and not overwrite_preprocessed:
            recording = si.load_extractor(preprocessed_path)
        else:
            recording = se.read_spikeglx(probe_path, stream_id=f'imec{probe_id}.ap')
            if config['artifact_correction']['do']:
                artifact_times = run_tprime_alignment(
                    input_dir=input_dir,
                    probe_id=probe_id,
                    tprime_config=config['tprime']
                )
                # recording = artifact_correction(
                #     recording=recording,
                #     artifact_times=artifact_times,
                #     window_ms=config['artifact_correction']['window_ms'],
                # )

                artifact_times_frame = artifact_times * recording.sampling_frequency
                recording = si.remove_artifacts(
                    recording, 
                    list_triggers=artifact_times_frame,
                    ms_before=config['artifact_correction']['ms_before'],
                    ms_after=config['artifact_correction']['ms_after'],
                    mode=config['artifact_correction']['mode'],
                )

            recording = recording.save(
                folder = preprocessed_path, 
                format='binary', 
                overwrite=True,
                **config['sorters']['job_kwargs'],
            )

        for sorter in sorter_configs['sorter_list']:
            folder = Path(probe_path) / sorter
            # if (folder / 'phy').exists():
            #     logger.info(f'Sorting already done for {sorter} on IMEC probe {probe_id}')
            #     continue

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

            ### The code below gives wrong output in phy, even with the output of Kilosort 
            ### looking fine, I dont know why. All clusters look like noise and 
            ### located at the tip of the probe

            # sorting = si.read_kilosort(folder)
            # # sorting = si.read_sorter_folder(folder)
            # sorting = si.remove_duplicated_spikes(
            #     sorting=sorting, 
            #     censored_period_ms=0.3,
            #     )
            # logger.info('Done running Kilosort for IMEC probe {}.'.format(probe_id))

            # sorting_analyzer = si.create_sorting_analyzer(
            #     sorting=sorting, 
            #     recording=recording, 
            #     sparse=True,
            #     overwrite=True,
            #     **config['sorters']['job_kwargs'],
            #     )

            # _ = sorting_analyzer.compute('random_spikes', method="uniform", max_spikes_per_unit=500)
            # _ = sorting_analyzer.compute('waveforms',  ms_before=1.5,ms_after=2., **config['sorters']['job_kwargs'])
            # _ = sorting_analyzer.compute('templates', operators=["average", "median", "std"])
            # _ = sorting_analyzer.compute('noise_levels')

            # # _ = sorting_analyzer.compute('spike_amplitudes', n_jobs = config['sorters']['job_kwargs']['n_jobs'])
            # # _ = sorting_analyzer.compute("spike_locations", n_jobs = config['sorters']['job_kwargs']['n_jobs'])

            # # _ = sorting_analyzer.compute(['correlograms', 'template_similarity', 'quality_metrics'],
            # #              extension_params=dict(quality_metrics=dict(metric_names=['snr', 'isi_violation', 'presence_ratio']))
            # #              )
            
            # # sorting_analyzer.save_as(
            # #     format='zarr', 
            # #     folder= folder / 'sorting_analyzer.zarr',
            # #     )
            
            # si.export_to_phy(
            #     sorting_analyzer=sorting_analyzer, 
            #     output_folder=folder / 'phy', 
            #     compute_pc_features=False,
            #     remove_if_exists=True,
            #     copy_binary=False,
            #     # use_relative_path=True,
            #     **config['sorters']['job_kwargs'],
            #     )
            # si.export_report(
            #     sorting_analyzer=sorting_analyzer, 
            #     output_folder=folder / 'report',
            #     remove_if_exists=True,
            #     **config['sorters']['job_kwargs'],
            #     )
            
        # shutil.rmtree(Path(probe_path) / 'preprocess', ignore_errors=True)

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