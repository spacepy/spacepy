:: Remove the Windows build system
@ECHO OFF

"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py37
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py38
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py39
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py310
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py311
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py312
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py313
start /wait "" "%SYSTEMDRIVE%\Miniconda3\Uninstall-Miniconda3.exe" /S /D=%SYSTEMDRIVE%\Miniconda3
