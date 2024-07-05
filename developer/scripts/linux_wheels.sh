#!/usr/bin/env bash

# ONLY for building local wheels, which are not tagged as manylinux compatible
# Get miniconda if you need it
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda
VERSIONS=(
    "3.6"
    "3.7"
    "3.8"
    "3.9"
    "3.10"
    "3.11"
    "3.12"
)
for PYVER in "${VERSIONS[@]}"
do
    ENVNAME=spacepy${PYVER//.}
    ~/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
    source ~/miniconda/bin/activate ${ENVNAME}
    if [ ${PYVER} = "3.6" ]; then
	conda install -y wheel
	pip install build
	BUILD=pyproject-build
    elif [ ${PYVER} = "3.9" ]; then
	conda install -y python-build wheel
	BUILD=python-build
    elif [ ${PYVER} = "3.12" ]; then
	conda install -y python-build wheel
	BUILD=python-build
    else
	conda install -y python-build wheel
	BUILD=python-build
    fi
    rm -rf build
    PYTHONNOUSERSITE=1 PYTHONPATH= SPACEPY_RELEASE=1 ${BUILD} -w -n -x .
    conda deactivate
    ~/miniconda/bin/conda env remove -y --name ${ENVNAME}
done
#rm -rf ~/miniconda
#rm ./Miniconda3-latest-Linux-x86_64.sh

