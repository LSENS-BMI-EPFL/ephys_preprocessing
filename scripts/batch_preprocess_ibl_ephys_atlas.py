#! /usr/bin/env python3
"""
@author: Axel Bisi
@project: ephys_utils
@file: batch_preprocess_ibl_ephys_atlas.py
@time: 8/19/2025
"""

import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import datetime

# ---------------- CONFIG ----------------
BASE_DIR = Path(r"M:\analysis\Axel_Bisi\data")
SCRIPT_PATH = Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_ibl_ephys_atlas.py")
CONFIG_FILE = Path(r"C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml")
CONDA_ENV = "iblenv"
MAX_WORKERS = 3   # max jobs in parallel if --parallel is used
LOG_DIR = Path("logs")
LOG_DIR.mkdir(exist_ok=True)

# List of mouse_id inputs you want to process
INPUTS = ["AB116", "AB120"]
# -----------------------------------------


def find_sessions(mouse_id: str):
    """Find session folders that contain an Ephys folder with catgt_*."""
    mouse_dir = BASE_DIR / mouse_id
    if not mouse_dir.exists():
        print(f"[WARNING] {mouse_id} not found under {BASE_DIR}")
        return []

    sessions = []
    for session in mouse_dir.iterdir():
        ephys_dir = session / "Ephys"
        if ephys_dir.exists():
            for catgt in ephys_dir.glob("catgt_*"):
                if catgt.is_dir():
                    sessions.append(catgt)
    return sessions


def run_job(mouse_id: str, catgt_path: Path, count: int):
    """Run one preprocessing job via conda run, logging output."""
    log_file = LOG_DIR / f"batchdriver_{mouse_id}_{count}.txt"

    # Command for conda run
    cmd = [
        #"conda", "run", "-n", CONDA_ENV,
        "python", str(SCRIPT_PATH),
        "--input", str(catgt_path),
        "--config", str(CONFIG_FILE)
    ]

    # Force UTF-8 logging for subprocess
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"

    # Force script to write its own logs inside LOG_DIR
    cwd = LOG_DIR

    with open(log_file, "w", encoding="utf-8", errors="replace") as lf:
        lf.write(f"[{datetime.datetime.now()}] Starting job for {mouse_id} at {catgt_path}\n")
        lf.flush()
        process = subprocess.Popen(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=cwd,   # 👈 ensures called script’s logs also go into LOG_DIR
            env=env
        )
        process.wait()
        lf.write(f"[{datetime.datetime.now()}] Finished job for {mouse_id} "
                 f"with return code {process.returncode}\n")

    return mouse_id, catgt_path, log_file, process.returncode


def main(parallel: bool):
    all_jobs = []
    for mouse_id in INPUTS:
        sessions = find_sessions(mouse_id)
        if not sessions:
            print(f"[WARNING] No sessions with catgt_* found for {mouse_id}")
            continue

        for i, session in enumerate(sessions):
            all_jobs.append((mouse_id, session, i))

    if not all_jobs:
        print("No jobs to run.")
        return

    if parallel:
        print(f"Launching {len(all_jobs)} jobs with up to {MAX_WORKERS} in parallel...")

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(run_job, mouse, path, i): (mouse, path)
                       for mouse, path, i in all_jobs}

            for future in as_completed(futures):
                mouse_id, path = futures[future]
                try:
                    result = future.result()
                    print(f"[OK] {result[0]} at {result[1]} → "
                          f"log: {result[2]} (exit {result[3]})")
                except Exception as e:
                    print(f"[ERROR] {mouse_id} at {path}: {e}")

    else:
        print(f"Launching {len(all_jobs)} jobs sequentially...")
        for mouse_id, path, i in all_jobs:
            try:
                result = run_job(mouse_id, path, i)
                print(f"[OK] {result[0]} at {result[1]} → "
                      f"log: {result[2]} (exit {result[3]})")
            except Exception as e:
                print(f"[ERROR] {mouse_id} at {path}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch preprocess ephys data")
    parser.add_argument("--parallel", action="store_true",
                        help="Run jobs in parallel (default: sequential)")
    args = parser.parse_args()

    main(parallel=args.parallel)
