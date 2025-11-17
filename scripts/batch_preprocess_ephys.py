#! /usr/bin/env python3
"""
Batch runner for IBL ephys preprocessing.

Features:
- Multi-script processing
- Per-script session finders
- loguru colored logs + master log file
- per-mouse/per-script log folder structure
- tqdm progress bars (per-mouse over scripts)
- CSV summary report with job details and runtimes
"""

import os
import sys
import yaml
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import datetime
import csv
from loguru import logger
from tqdm import tqdm
import platform

# ---------------- CONFIG ----------------

EXPERIMENTER = 'Axel_Bisi'
if EXPERIMENTER == 'Axel_Bisi':
    machine = platform.node()
    if machine == 'SV-07-014':
        CONFIG_FILE = r'C:\Users\bisi\Github\ephys_preprocessing\preprocessing\preprocess_config.yaml'
        SCRIPTS = {
            "preprocess_spikesort": Path(r"C:\Users\bisi\Github\ephys_preprocessing\preprocessing\preprocess_spikesort.py"),
            "preprocess_sync": Path(r"C:\Users\bisi\Github\ephys_preprocessing\preprocessing\preprocess_sync.py"),
            # "preprocess_ibl_ephys_atlas": Path(...),
        }
    elif machine == 'SV-07-081':
        CONFIG_FILE = r'C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml'
        SCRIPTS = {
            "preprocess_spikesort": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_spikesort.py"),
            "preprocess_sync": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_sync.py"),
            # "preprocess_ibl_ephys_atlas": Path(...),
        }
    else:
        # default in this branch
        CONFIG_FILE = r'C:\Users\bisi\Github\ephys_preprocessing\preprocessing\preprocess_config.yaml'
        SCRIPTS = {
            "preprocess_spikesort": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_spikesort.py"),
            "preprocess_sync": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_sync.py"),
        }
else:
    CONFIG_FILE = Path(r"C:\Users\bisi\Github\ephys_preprocessing\preprocessing\preprocess_config.yaml")
    SCRIPTS = {
        "preprocess_spikesort": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_spikesort.py"),
        "preprocess_sync": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_sync.py"),
        "preprocess_ibl_ephys_atlas": Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_ibl_ephys_atlas.py"),
    }

# Load config file as dict
with open(CONFIG_FILE, 'r', encoding='utf8') as stream:
    CONFIG = yaml.safe_load(stream)

RAW_DIR = Path(CONFIG['raw_data_path'])
BASE_DIR = Path(CONFIG['output_path'])

LOG_ROOT = Path("logs")
LOG_ROOT.mkdir(exist_ok=True)

MAX_WORKERS = 6

# Mice to process
INPUTS = ["AB150", "AB151", "AB152", "AB153", "AB154", "AB155", "AB156", "AB157", "AB158", "AB159",
          "AB160", "AB161", "AB162", "AB163", "AB164"]

# For testing override a subset
# INPUTS = ["AB150"]
# -----------------------------------------

# LOGURU SETUP
logger.remove()
logger.add(lambda msg: print(msg, end=""), colorize=True,
           format="<green>{time:HH:mm:ss}</green> | <level>{level}</level> | {message}")
logger.add(LOG_ROOT / "master_log.txt", rotation="5 MB", backtrace=True, diagnose=False)


# -------- JOB DISCOVERY HELPERS --------

def find_processed_sessions(mouse_id: str):
    """
    Find processed session folders for a mouse.
    Returns catgt_* paths inside each session's Ephys folder.
    Example returned path:
      M:\analysis\Axel_Bisi\data\<mouse_id>\<session_id>\Ephys\catgt_...
    """
    mouse_dir = BASE_DIR / mouse_id
    if not mouse_dir.exists():
        logger.warning(f"{mouse_id} not found under {BASE_DIR}")
        return []

    sessions = []
    for session in mouse_dir.iterdir():
        ephys_dir = session / "Ephys"
        if ephys_dir.exists() and ephys_dir.is_dir():
            for catgt in ephys_dir.glob("catgt_*"):
                if catgt.is_dir():
                    sessions.append(catgt)
    return sessions


def find_raw_sessions(mouse_id: str):
    """
    Find raw session folders for a mouse.
    Returns Ephys folders under RAW_DIR/<mouse_id>/Recording/<session>
    Example:
      M:\data\<mouse_id>\Recording\<session_id>\Ephys
    """
    mouse_dir = RAW_DIR / mouse_id / "Recording"
    if not mouse_dir.exists() or not mouse_dir.is_dir():
        logger.warning(f"{mouse_id}/Recording not found under {RAW_DIR}")
        return []

    sessions = []
    for session in mouse_dir.iterdir():
        ephys_dir = session / "Ephys"
        if ephys_dir.exists() and ephys_dir.is_dir():
            sessions.append(ephys_dir)
    return sessions


# You can add more finders here if needed


# -------- JOB EXECUTION --------

def run_job(mouse_id: str, session_path: Path, count: int, script_key: str, utils_path: Path):
    """Run preprocessing job for one script + one session and return result dict."""

    script_path = SCRIPTS[script_key]

    # Validate script file exists
    if not script_path.exists() or not script_path.is_file():
        logger.error(f"Script for key '{script_key}' not found at {script_path}")
        return {
            "timestamp": datetime.datetime.now().isoformat(),
            "mouse_id": mouse_id,
            "script": script_key,
            "session_path": str(session_path),
            "log_file": "",
            "return_code": -1,
            "status": "error",
            "runtime_sec": 0.0,
            "error": f"Script not found: {script_path}"
        }

    # Per-mouse/per-script log folder
    log_dir = LOG_ROOT / mouse_id / script_key
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{mouse_id}_{script_key}_{count}.txt"

    # Use the same Python interpreter as the runner
    python_exec = sys.executable

    cmd = [
        python_exec,
        str(script_path),
        "--input", str(session_path),
        "--config", str(CONFIG_FILE),
    ]

    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"

    # Add utils folder to PYTHONPATH for subprocess
    if utils_path:
        utils_path = Path(utils_path).resolve()
        if "PYTHONPATH" in env and env["PYTHONPATH"]:
            env["PYTHONPATH"] = str(utils_path) + os.pathsep + env["PYTHONPATH"]
        else:
            env["PYTHONPATH"] = str(utils_path)

    # Run the script in its folder (so relative imports in the script work)
    cwd = script_path.parent

    start_time = datetime.datetime.now()

    with open(log_file, "w", encoding="utf-8", errors="replace") as lf:
        lf.write(f"[{start_time}] Starting job for {mouse_id}, script={script_key}, path={session_path}\n")
        lf.flush()

        process = subprocess.Popen(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=str(cwd),
            env=env
        )
        process.wait()

        end_time = datetime.datetime.now()
        lf.write(f"[{end_time}] Finished job for {mouse_id} (script={script_key}) with return code {process.returncode}\n")

    runtime = (end_time - start_time).total_seconds()

    return {
        "timestamp": end_time.isoformat(),
        "mouse_id": mouse_id,
        "script": script_key,
        "session_path": str(session_path),
        "log_file": str(log_file),
        "return_code": process.returncode,
        "status": "success" if process.returncode == 0 else "error",
        "runtime_sec": runtime,
    }


# ---------------- PROCESS ONE MOUSE ----------------

def process_mouse(mouse_id: str, scripts_sessions: dict, utils_path: Path = None):
    """
    Run scripts sequentially for one mouse.
    scripts_sessions: dict mapping script_key -> list of session Paths
    Returns list of result dicts for every (script, session) executed.
    """

    logger.info(f"[{mouse_id}] Processing started")
    results = []

    # iterate scripts that have sessions (preserve SCRIPTS order)
    scripts_to_run = [k for k in SCRIPTS.keys() if k in scripts_sessions and scripts_sessions[k]]

    # progress bar over scripts for this mouse
    for script_key in tqdm(scripts_to_run, desc=f"[{mouse_id}] Scripts", unit="script", leave=False):
        sessions = scripts_sessions.get(script_key, [])
        logger.info(f"[{mouse_id}] Running script: {script_key} on {len(sessions)} session(s)")

        for i, session in enumerate(sessions):
            logger.info(f"[{mouse_id}]   Session {i+1}/{len(sessions)}: {session}")
            res = run_job(mouse_id, session, i, script_key, utils_path)
            results.append(res)

            if res["return_code"] != 0:
                logger.error(f"[{mouse_id}] ERROR in script={script_key}, session={session} (exit {res['return_code']})")
                # continue to next session/script (no abort) — change behavior here if you want to stop on error

    logger.success(f"[{mouse_id}] All scripts completed")
    return results


# ---------------- MAIN PIPELINE ----------------

def main(parallel: bool):
    # SESSION_FINDERS mapping: script_key -> finder function
    SESSION_FINDERS = {
        "preprocess_spikesort": find_raw_sessions,
        "preprocess_sync": find_processed_sessions,
        "preprocess_ibl_ephys_atlas": find_processed_sessions
    }

    # Path to utils (adjust if different)
    default_utils_path = Path(__file__).parent.parent / "utils"
    if not default_utils_path.exists():
        # fallback to likely repo location used earlier
        default_utils_path = Path(r"C:\Users\bisi\Github\ephys_preprocessing\utils")
    default_utils_path = default_utils_path.resolve()

    # Build one mouse job: for each mouse, build a dict script->sessions
    mouse_jobs = []
    for mouse_id in INPUTS:
        scripts_sessions = {}
        any_sessions_found = False

        for script_key in SCRIPTS.keys():
            finder = SESSION_FINDERS.get(script_key, find_processed_sessions)
            try:
                sessions = finder(mouse_id) or []
            except Exception as e:
                logger.error(f"Finder for {script_key} raised error for {mouse_id}: {e}")
                sessions = []

            if sessions:
                any_sessions_found = True
                scripts_sessions[script_key] = sessions

        if not any_sessions_found:
            logger.warning(f"No sessions found for {mouse_id} for any script; skipping mouse.")
            continue

        mouse_jobs.append((mouse_id, scripts_sessions))

    if not mouse_jobs:
        logger.info("No jobs to run.")
        return

    logger.info('Sessions to process (summary):')
    for mouse_id, scripts_sessions in mouse_jobs:
        for sk, sess_list in scripts_sessions.items():
            logger.info(f" - Mouse: {mouse_id}, Script: {sk}, Sessions: {len(sess_list)}")

    results = []

    # ---- RUN MICE IN PARALLEL ----
    if parallel:
        logger.info(f"Running {len(mouse_jobs)} mice in parallel with {MAX_WORKERS} workers...")

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {}
            for mouse_id, scripts_sessions in mouse_jobs:
                fut = executor.submit(process_mouse, mouse_id, scripts_sessions, default_utils_path)
                futures[fut] = mouse_id

            for fut in tqdm(as_completed(futures), total=len(futures), desc="Mice", unit="mouse"):
                mouse_id = futures[fut]
                try:
                    mouse_results = fut.result()
                    results.extend(mouse_results)
                    logger.success(f"[{mouse_id}] finished (added {len(mouse_results)} job results)")
                except Exception as e:
                    logger.error(f"[{mouse_id}] FAILURE: {e}")

    # ---- RUN MICE SEQUENTIALLY ----
    else:
        logger.info("Running mice sequentially...")
        for mouse_id, scripts_sessions in tqdm(mouse_jobs, desc="Mice", unit="mouse"):
            try:
                mouse_results = process_mouse(mouse_id, scripts_sessions, default_utils_path)
                results.extend(mouse_results)
                logger.success(f"[{mouse_id}] finished (added {len(mouse_results)} job results)")
            except Exception as e:
                logger.error(f"[{mouse_id}] FAILURE: {e}")

    # ---- CSV SUMMARY ----
    if results:
        csv_path = LOG_ROOT / "batch_summary.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        logger.success(f"CSV summary written → {csv_path}")
    else:
        logger.info("No results to write to CSV.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch preprocess ephys data"
    )
    parser.add_argument("--parallel", action="store_true",
                        help="Run mice in parallel (default: sequential)")
    args = parser.parse_args()

    main(parallel=args.parallel)
