:: Remove the Windows build system
@ECHO OFF

%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py27_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py27_64
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py36_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py36_64
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py37_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py37_64
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py38_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py38_64
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py39_32
%USERPROFILE%\Miniconda3\Scripts\conda env remove -y --name py39_64
start /wait "" %USERPROFILE%\Miniconda3\Uninstall-Miniconda3.exe /S /D=%USERPROFILE%\Miniconda3
