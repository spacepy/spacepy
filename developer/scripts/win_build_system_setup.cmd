:: Set up a Windows build system
:: Download 64-bit Python 3.7 miniconda from
:: https://docs.conda.io/en/latest/miniconda.html
:: Put the updates in here, too
@ECHO OFF

set CONDA_PKGS_DIRS=%USERPROFILE%\Miniconda3\PKGS64
set CONDA_SUBDIR=win-64
set CONDA_FORCE_32_BIT=

start /wait "" %USERPROFILE%\Downloads\Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /NoRegistry=1 /S /D=%USERPROFILE%\Miniconda3
:: base environment needs to be activated
CALL %USERPROFILE%\Miniconda3\Scripts\activate
CALL conda update -y conda
CALL conda create -y -n py3_64 python=3.6
CALL conda create -y -n py2_64 python=2.7

set CONDA_PKGS_DIRS=%USERPROFILE%\Miniconda3\PKGS32
set CONDA_SUBDIR=win-32
set CONDA_FORCE_32_BIT=1
CALL conda create -y -n py3_32 python=3.6
CALL conda create -y -n py2_32 python=2.7
CALL conda deactivate

FOR %%B in (32 64) DO (FOR %%P in (2 3) DO CALL :installs %%B %%P)

GOTO :EOF

:installs
IF "%1"=="32" (
    set CONDA_PKGS_DIRS=%USERPROFILE%\Miniconda3\PKGS32
    set CONDA_SUBDIR=win-32
    set CONDA_FORCE_32_BIT=1
) ELSE (
    set CONDA_PKGS_DIRS=%USERPROFILE%\Miniconda3\PKGS64
    set CONDA_SUBDIR=win-64
    set CONDA_FORCE_32_BIT=
)
CALL %USERPROFILE%\Miniconda3\Scripts\activate py%2_%1
CALL conda install -y numpy scipy matplotlib networkx m2w64-gcc-fortran libpython h5py
:: libpython sets things up to use ming by default, otherwise try distutils.cfg
:: note we need libpython or else ffnet requires MSVC to install
SET FC_VENDOR=gfortran
CALL pip install ffnet
CALL conda deactivate
::This turns off echo!
GOTO :EOF
