:: Set up a Windows build system
:: Download 64-bit Python 3.7 miniconda from
:: https://docs.conda.io/en/latest/miniconda.html
:: Put the updates in here, too
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
CALL conda create -y -n py39_64 python=3.9
CALL conda create -y -n py38_64 python=3.8
CALL conda create -y -n py37_64 python=3.7
CALL conda create -y -n py36_64 python=3.6
CALL conda create -y -n py27_64 python=2.7

set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS32
set CONDA_SUBDIR=win-32
set CONDA_FORCE_32_BIT=1
CALL conda create -y -n py39_32 python=3.9
CALL conda create -y -n py38_32 python=3.8
CALL conda create -y -n py37_32 python=3.7
CALL conda create -y -n py36_32 python=3.6
CALL conda create -y -n py27_32 python=2.7
CALL conda deactivate

IF "%1"=="build" (
    set ACTION=BUILD
) ELSE (
    set ACTION=TEST
)

FOR %%B in (32 64) DO (FOR %%P in (27 36 37 38 39) DO CALL :installs %%B %%P)

GOTO :EOF

:installs
IF "%1"=="32" (
    set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS32
    set CONDA_SUBDIR=win-32
    set CONDA_FORCE_32_BIT=1
) ELSE (
    set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS64
    set CONDA_SUBDIR=win-64
    set CONDA_FORCE_32_BIT=
)
CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%2_%1

:: Are we building SpacePy wheels, or testing?
IF "%ACTION%"=="BUILD" (
    :: Get the compiler
    CALL conda install -y m2w64-gcc-fortran libpython
    set NUMPY="numpy"
    :: minimum version for each Python version
    IF "%2"=="39" (
    :: 1.18 works on 3.9, but there's no Windows binary wheel.
    :: 1.19.4 has Win10 2004 bug on 64-bit, but
    :: we'll avoid it on 32-bit as well, no point getting picky...
        set NUMPY="numpy>=1.19.5,<1.20.0"
    )
    IF "%2"=="38" (
        set NUMPY="numpy>=1.17.0,<1.19.0"
    )
    IF "%2"=="37" (
        set NUMPY="numpy>=1.15.1,<1.16.0"
    )
    IF "%2"=="36" (
    :: 1.12 works on 3.6, but f2py fails on Windows
        set NUMPY="numpy>=1.13.0,<1.14.0"
    )
    IF "%2"=="35" (
        set NUMPY="numpy>=1.10.0,<1.11.0"
    )
    IF "%2"=="27" (
        set NUMPY="numpy>=1.10.0,<1.11.0"
    )
    CALL pip install !NUMPY!
) ELSE (
    :: Testing. Get the latest of everything
    IF "%2"=="39" (
        CALL pip install numpy
        CALL pip install numpy scipy matplotlib networkx h5py ffnet astropy
    ) ELSE (
        CALL conda install -y numpy scipy matplotlib networkx h5py astropy
    )
    :: libpython sets things up to use ming by default, otherwise try distutils.cfg
    :: note we need libpython or else ffnet requires MSVC to install
    SET FC_VENDOR=gfortran
    CALL pip install ffnet
)
CALL conda deactivate
::This turns off echo!
GOTO :EOF
