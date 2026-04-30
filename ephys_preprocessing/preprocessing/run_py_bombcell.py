import os
import pathlib
from loguru import logger
import bombcell as bc
import yaml

from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording, extract_ks_version
from ephys_preprocessing.utils.phylib_utils import load_phy_model
def main(input_dir, config):
    input_dir = os.path.join(input_dir, [f for f in os.listdir(input_dir) if 'catgt' in f][0])
    catgt_epoch_name = os.path.basename(input_dir)
    epoch_name = catgt_epoch_name.lstrip('catgt_')

    probe_folders = [f for f in os.listdir(input_dir) if 'imec' in f]
    probe_ids = [f[-1] for f in probe_folders]

    # Perform computations for each probe separately
    for probe_id in probe_ids:

        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[0]
        if not check_if_valid_recording(config, mouse_id, probe_id):
            continue

        probe_folder = '{}_imec{}'.format(epoch_name, probe_id)
        probe_path = os.path.join(input_dir, probe_folder)

        kilosort_folders = pathlib.Path(probe_path).glob('kilosort*')
        for kilosort_folder in kilosort_folders:
            # Set paths
            ks_name = kilosort_folder.name
            kilosort_version = extract_ks_version(ks_name)
            if mouse_id.startswith('AB') or mouse_id.startswith('MH'):
                kilosort_path = os.path.join(kilosort_folder)
            else:
                kilosort_path = os.path.join(kilosort_folder, 'sorter_output')

            save_path = os.path.join(kilosort_path, 'bombcell')

            try:
                apbin_fname = '{}_tcat_corrected.imec{}.ap.bin'.format(epoch_name, probe_id)
                meta_fname = '{}_tcat_corrected.imec{}.ap.meta'.format(epoch_name, probe_id)
                path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)
                path_to_meta = os.path.join(input_dir, probe_folder, meta_fname)

            except FileNotFoundError:
                apbin_fname = '{}_tcat.imec{}.ap.bin'.format(epoch_name, probe_id)
                meta_fname = '{}_tcat.imec{}.ap.meta'.format(epoch_name, probe_id)
                path_to_apbin = os.path.join(input_dir, probe_folder, apbin_fname)
                path_to_meta = os.path.join(input_dir, probe_folder, meta_fname)

            param = bc.get_default_parameters(
                kilosort_path, 
                raw_file=path_to_apbin,
                meta_file=path_to_meta,
                kilosort_version=kilosort_version,
            )

            logger.info('Running bombcell for IMEC probe {}.'.format(probe_id))

            # Compute quality metrics
            quality_metrics, param, unit_type, unit_type_string, = bc.run_bombcell(
                kilosort_path, save_path, param,
            )

            # Compute ephys properties for cell type classification
            ephys_param = bc.get_ephys_parameters(kilosort_path)

            # Compute all ephys properties - now defaults to ks_dir/bombcell
            ephys_properties, ephys_param = bc.run_all_ephys_properties(kilosort_path, ephys_param, save_path=save_path)

        # cluster_info table creation
        logger.info('Creating cluster_info table for IMEC probe {}.'.format(probe_id))
        phy_model = load_phy_model(os.path.join(kilosort_path, 'params.py'))
        phy_model.create_metrics_dataframe()
        phy_model.save_metrics_tsv(os.path.join(kilosort_path, 'cluster_info.tsv'))

    return

if __name__ == '__main__':
    input_dir = '/home/lebert/lsens_srv/analysis/Jules_Lebert/data/JL007/JL007_20250603_150143/Ephys/'
    config_file = '/home/lebert/code/spikesorting_pipeline/spikeinterface_preprocessing/ephys_preprocessing/scripts/preprocess_config_si.yaml'
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    main(input_dir, config)