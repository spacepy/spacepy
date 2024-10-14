#!/usr/bin/env zsh

# Get miniconda if you need it
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
#bash ./Miniconda3-latest-MacOSX-arm64.sh -b -p ~/opt/miniconda
PYVER=3.10
ENVNAME=spacepy${PYVER//.}
~/opt/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
source ~/opt/miniconda/bin/activate ${ENVNAME}
conda install -y python-build wheel gfortran~=11.2.0
BUILD=python-build
rm -rf build
PYTHONNOUSERSITE=1 PYTHONPATH= SPACEPY_RELEASE=1 SDKROOT=/opt/MacOSX11.0.sdk ${BUILD} -w -n -x .
conda deactivate
# if these vars aren't cleared, next conda env will use them for compilers
badvars=$(env | sed -n 's/^\(CONDA_BACKUP_[^=]*\)=.*/\1/p')
unset $(echo $badvars)
export $(echo $badvars)
~/opt/miniconda/bin/conda env remove -y --name ${ENVNAME}
#rm -rf ~/opt/miniconda
#rm ./Miniconda3-latest-MacOSX-arm64.sh
