:: Set up a Windows build system
:: Download 64-bit Python 3.7 miniconda from
:: https://docs.conda.io/en/latest/miniconda.html
:: Put the updates in here, too
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

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

IF "%1"=="build" (
    set ACTION=BUILD
) ELSE (
    set ACTION=TEST
)

FOR %%P in (37 38 39 310 311 312) DO CALL :installs %%P

GOTO :EOF

:installs

CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%1

:: Are we building SpacePy wheels, or testing?
IF "%ACTION%"=="BUILD" (
    :: Get the compiler
    CALL conda install -y m2w64-gcc-fortran libpython
    CALL conda install -y python-build wheel
) ELSE (
    :: Testing. Get the latest of everything
    CALL conda install -y numpy scipy matplotlib h5py astropy
)
CALL conda deactivate
::This turns off echo!
GOTO :EOF
