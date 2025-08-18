@echo on
REM ==========================================
REM batch_preprocess_parallel_test.bat
REM Verbose test run: process only the first 3 "catgt*" folders
REM ==========================================

REM --- CONFIGURATION ---
set "BASE_DIR=M:\analysis\Axel_Bisi\data"
set "SCRIPT_PATH=C:\Users\bisi\ephys_utils\preprocess_ibl_ephys_atlas.py"
set "CONFIG_FILE=C:\Users\bisi\ephys_utils\preprocessing\preprocess_config.yaml"
set "CONDA_ENV=iblenv"

REM List of AB folders you want to process
set INPUTS=AB123 AB164 AB163

echo Base directory  : %BASE_DIR%
echo Script path     : %SCRIPT_PATH%
echo Config file     : %CONFIG_FILE%
echo Conda env       : %CONDA_ENV%
echo Inputs          : %INPUTS%
echo.

REM --- Loop over each AB input ---
setlocal enabledelayedexpansion
set COUNT=0

for %%A in (%INPUTS%) do (
    echo -----------------------------------------
    echo Processing input %%A
    echo -----------------------------------------

    REM --- Find the session folder containing Ephys ---
    set "SESSION_FOUND="
    for /D %%S in ("%BASE_DIR%\%%A\*") do (
        if exist "%%S\Ephys" (
            set "SESSION_FOUND=%%S"
            goto :FOUND_SESSION
        )
    )
    :FOUND_SESSION

    if defined SESSION_FOUND (
        echo Found session folder: !SESSION_FOUND!

        REM Find catgt_* folders inside Ephys
        for /D %%C in ("!SESSION_FOUND!\Ephys\catgt_*") do (
            set "INPUT_DIR=%%C"
            echo Using input dir: !INPUT_DIR!

            REM Launch Python script in parallel using conda run
            START "job!COUNT!" cmd /c ^
            ""C:\Users\bisi\AppData\Local\anaconda3\condabin\conda.bat" run -n %CONDA_ENV% python "%SCRIPT_PATH%" --input "!INPUT_DIR!" --config "%CONFIG_FILE%""

            set /a COUNT+=1
        )
    ) else (
        echo [WARNING] No session folder with Ephys found for %%A
    )
)

echo.
echo Launched !COUNT! jobs in parallel.
echo.
pause