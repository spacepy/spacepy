:: Remove the Windows build system
@ECHO OFF

"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py27_32
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py27_64
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py36_32
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py36_64
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py37_32
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py37_64
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py38_32
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py38_64
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py39_32
"%SYSTEMDRIVE%\Miniconda3\Scripts\conda" env remove -y --name py39_64
start /wait "" "%SYSTEMDRIVE%\Miniconda3\Uninstall-Miniconda3.exe" /S /D=%SYSTEMDRIVE%\Miniconda3
