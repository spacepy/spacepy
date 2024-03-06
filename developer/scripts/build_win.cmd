:: Build spacepy on windows
:: Assume win_build_system_setup.cmd has already been run
@ECHO OFF
SETLOCAL EnableDelayedExpansion

:: Get any system Python out of the way
set PYTHONPATH=
set PATH=

FOR %%P in (36 37 38 39 310 311 312) DO CALL :build %%P

GOTO :EOF

:build

CALL "%SYSTEMDRIVE%\Miniconda3\Scripts\activate" py%1
pushd %~dp0\..\..\
rmdir /s /q build 2> nul
IF "%1"=="36" (
    CALL pyproject-build -w -n -x
) ELSE (
    CALL python-build -w -n -x
)
popd
::This turns off echo!
CALL conda deactivate
GOTO :EOF
