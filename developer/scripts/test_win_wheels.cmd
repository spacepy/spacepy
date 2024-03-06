:: Run unit tests on all windows wheels, run from unittest dir
:: Download 64-bit miniconda from
:: https://docs.conda.io/en/latest/miniconda.html
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS64
set CONDA_SUBDIR=win-64
set CONDA_FORCE_32_BIT=

start /wait "" "%USERPROFILE%\Downloads\Miniconda3-latest-Windows-x86_64.exe" /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /NoRegistry=1 /S /D=%SYSTEMDRIVE%\Miniconda3
:: base environment needs to be activated
CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate"
CALL conda update -y conda
CALL conda create -y -n py312 python=3.12
CALL conda create -y -n py311 python=3.11
CALL conda create -y -n py310 python=3.10
CALL conda create -y -n py39 python=3.9
CALL conda create -y -n py38 python=3.8
CALL conda create -y -n py37 python=3.7
CALL conda create -y -n py36 python=3.6


FOR %%P in (36 37 38 39 310 311 312) DO CALL :installs %%P

GOTO :EOF

:installs
set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS64
set CONDA_SUBDIR=win-64
set CONDA_FORCE_32_BIT=

CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%1
::CALL pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ --only-binary spacepy "spacepy>0.0d0"
:: Installs latest prerelease while keeping deps stable
CALL pip install --find-links ..\dist --only-binary spacepy "spacepy>0.0d0"
CALL python test_all.py > test_output_py%1.txt 2>&1
CALL conda deactivate
::This turns off echo!
GOTO :EOF
