:: Build spacepy on windows
:: Assume win_build_system_setup.cmd has already been run
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

FOR %%P in (310) DO CALL :build %%P

GOTO :EOF

:build

CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%1
pushd %~dp0\..\..\
rmdir /s /q build 2> nul
CALL python-build -w -n -x

popd
::This turns off echo!
CALL conda deactivate
GOTO :EOF
