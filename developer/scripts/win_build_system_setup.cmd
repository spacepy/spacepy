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
CALL conda create -y -n py312_64 python=3.12
CALL conda create -y -n py311_64 python=3.11
CALL conda create -y -n py310_64 python=3.10
CALL conda create -y -n py39_64 python=3.9
CALL conda create -y -n py38_64 python=3.8
CALL conda create -y -n py37_64 python=3.7
CALL conda create -y -n py36_64 python=3.6

IF "%1"=="build" (
    set ACTION=BUILD
) ELSE (
    set ACTION=TEST
)

FOR %%B in (64) DO (FOR %%P in (36 37 38 39 310) DO CALL :installs %%B %%P)

GOTO :EOF

:installs
set CONDA_PKGS_DIRS=%SYSTEMDRIVE%\Miniconda3\PKGS64
set CONDA_SUBDIR=win-64
set CONDA_FORCE_32_BIT=

CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%2_%1

:: Are we building SpacePy wheels, or testing?
IF "%ACTION%"=="BUILD" (
    :: Get the compiler
    CALL conda install -y m2w64-gcc-fortran libpython
    set NUMPY="numpy"
    :: Build with the minimum version for each Python version
    IF "%2"=="312" (
        set NUMPY="numpy~=1.24.0"
    )
    IF "%2"=="311" (
        set NUMPY="numpy~=1.22.0"
    )
    IF "%2"=="310" (
        set NUMPY="numpy~=1.21.0"
    )
    IF "%2"=="39" (
    :: 1.18 works on 3.9, but there's no Windows binary wheel.
    :: 1.19.4 has Win10 2004 bug on 64-bit, but
    :: we'll avoid it on 32-bit as well, no point getting picky...
        set NUMPY="numpy~=1.19.5"
    )
    IF "%2"=="38" (
        set NUMPY="numpy~=1.17.0"
    )
    IF "%2"=="37" (
        set NUMPY="numpy~=1.15.1"
    )
    IF "%2"=="36" (
        set NUMPY="numpy~=1.15.1"
    )
    CALL pip install !NUMPY! build wheel
) ELSE (
    :: Testing. Get the latest of everything
    CALL conda install -y numpy scipy matplotlib h5py astropy
)
CALL conda deactivate
::This turns off echo!
GOTO :EOF
