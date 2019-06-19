:: Remove the Windows build system
@ECHO OFF

%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py2_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py2_64
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py3_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py3_64
start /wait "" %USERPROFILE%\Miniconda3\Uninstall-Miniconda3.exe /S /D=%USERPROFILE%\Miniconda3
