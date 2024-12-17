import os
from pathlib import Path
from loguru import logger
import yaml
import spikeinterface.extractors as se
import spikeinterface.full as si

from ephys_preprocessing.utils.ephys_utils import check_if_valid_recording

def main(input_dir, config):
    # epoch_name = os.listdir(input_dir)[0]
    epoch_name = list(input_dir.glob('catgt*'))[0].name
    probe_folders = list((input_dir / epoch_name).glob("*imec*"))
    # probe_folders = [f for f in os.listdir(os.path.join(input_dir, epoch_name)) if 'imec' in f]
    logger.info('Data to spike-sort: {}'.format(probe_folders))
    n_probes = len(probe_folders)

    for probe_id in range(n_probes):
        # Check if probe recording is valid
        mouse_id = epoch_name.split('_')[1]
        # if not check_if_valid_recording(config, mouse_id, probe_id):
        #     continue

        probe_folder = '{}_imec{}'.format(epoch_name.replace('catgt_', ''), probe_id)
        probe_path = os.path.join(input_dir, epoch_name, probe_folder)

        recording = se.read_spikeglx()

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