import sys
import os
import subprocess
import numpy as np

from pathlib import Path
from typing import Union, Optional
from loguru import logger

def run_tprime_alignment(
    input_dir: Union[str, Path],
    probe_id: int,
    tprime_config: dict,
    run_name: Optional[str] = None,
    force_rerun: bool = False,
) -> np.ndarray:
    """
    Run TPrime to get artifact times aligned to probe timebase.
    
    Parameters
    ----------
    input_dir : str or Path
        Directory containing CatGT-processed data
    probe_id : int
        Probe ID number
    tprime_config : dict
        Configuration dictionary containing:
        - tprime_path: path to TPrime executable
        - syncperiod: sync period for TPrime
    run_name : str or None
        Recording run name. If None, extracts from directory structure
    force_rerun : bool, default: False
        If True, runs TPrime even if output file exists. If False, loads existing output if available.
        
    Returns
    -------
    artifact_times : np.ndarray
        Array of artifact times in seconds aligned to probe timebase
    """
    input_dir = Path(input_dir)
    
    # Get run name if not provided
    if run_name is None:
        epoch_name = list(input_dir.glob("catgt*"))[0].name
        run_name = epoch_name[6:]
    
    # Find probe folder
    epoch_dir = input_dir / f"catgt_{run_name}"
    probe_folders = list(epoch_dir.glob(f"*imec{probe_id}*"))
    if not probe_folders:
        raise FileNotFoundError(f"No folder found for probe {probe_id}")
    probe_path = probe_folders[0]
    
    # Check for existing output file
    output_file = probe_path / f"whisker_stim_times_to_imec{probe_id}.txt"
    if output_file.exists() and not force_rerun:
        logger.info(f'Loading existing TPrime output for probe {probe_id}')
        return np.loadtxt(output_file)
    
    # Set up TPrime path based on OS
    if sys.platform.startswith('win'):
        tprime_exe = 'TPrime'
        shell = True
    elif sys.platform.startswith('linux'):
        tprime_exe = os.path.join(tprime_config['tprime_path'], "runit.sh").replace('\\', '/')
        shell = False
    else:
        raise NotImplementedError('OS not supported')
    
    # Get number of channels from meta file
    meta_files = list(probe_path.glob("*.ap.meta"))
    if not meta_files:
        raise FileNotFoundError("No .ap.meta file found")
    
    # Read number of channels from meta file
    with open(meta_files[0]) as f:
        for line in f:
            if line.startswith('nSavedChans'):
                n_channels = int(line.split('=')[1])
                break
    
    tostream_probe_edges_file = f'{run_name}_tcat.imec{probe_id}.ap.xd_{n_channels-1}_6_500.txt'
    nidq_stream_idx = 10 # arbitrary index number

    # Build and run TPrime command
    command = [
        tprime_exe,
        f'-syncperiod={tprime_config["syncperiod"]}',
        f'-tostream={probe_path / tostream_probe_edges_file}',
        f'-fromstream={nidq_stream_idx},{epoch_dir / f"{run_name}_tcat.nidq.xa_0_0.txt"}',
        f'-events={nidq_stream_idx},'
        f'{epoch_dir / f"{run_name}_tcat.nidq.xa_3_0.txt"},'
        f'{output_file}'
    ]
    
    logger.info(f'Running TPrime to sync whisker artifact times to IMEC probe {probe_id} timebase')
    subprocess.run(command, shell=shell, cwd=tprime_config['tprime_path'])
    
    # Read and return aligned artifact times
    if not output_file.exists():
        raise FileNotFoundError(f"TPrime failed to create output file: {output_file}")
        
    return np.loadtxt(output_file)