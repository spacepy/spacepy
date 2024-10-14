:: Build spacepy on windows
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

start /wait "" "%USERPROFILE%\Downloads\Miniconda3-latest-Windows-x86_64.exe" /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /NoRegistry=1 /S /D=%SYSTEMDRIVE%\Miniconda3
:: base environment needs to be activated
CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate"
CALL conda update -y conda
CALL conda create -y -n py310 python=3.10
CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py310
CALL conda install -y m2w64-gcc-fortran libpython
CALL conda install -y python-build wheel
pushd %~dp0\..\..\
rmdir /s /q build 2> nul
CALL python-build -w -n -x

popd
::This turns off echo!
CALL conda deactivate
