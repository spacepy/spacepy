:: Build spacepy on windows
:: Assume win_build_system_setup.cmd has already been run
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

FOR %%B in (64) DO (FOR %%P in (36 37 38 39 310) DO CALL :build %%B %%P)

GOTO :EOF

:build
set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS64
set CONDA_SUBDIR=win-64
set CONDA_FORCE_32_BIT=

set PYVER="%2"
CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%2_%1
pushd %~dp0\..\..\
rmdir /s /q build 2> nul
CALL pyproject-build -w -n -x
popd
::This turns off echo!
CALL conda deactivate
GOTO :EOF
