#!/usr/bin/env bash

# Get miniconda if you need it
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda
VERSIONS=(
    "3.6|numpy~=1.12.0"
    "3.7|numpy~=1.15.1"
    "3.8|numpy~=1.17.0"
    "3.9|numpy~=1.17.0"
    "3.10|numpy~=1.21.0"
    "3.11|numpy~=1.22.0"
    "3.12|numpy~=1.24.0"
)
for thisBuild in "${VERSIONS[@]}"
do
    IFS='|' read -ra tmpVar <<< "$thisBuild"
    PYVER=${tmpVar[0]}
    NUMPY="${tmpVar[1]}"
    ENVNAME=spacepy${PYVER//.}
    ~/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
    source ~/miniconda/bin/activate ${ENVNAME}
    if [ ${PYVER} = "3.6" ]; then
	conda install -y wheel ${NUMPY}
	pip install build
	BUILD=pyproject-build
    elif [ ${PYVER} = "3.9" ]; then
	conda install -y python-build wheel "cython<3"
	pip install --no-build-isolation ${NUMPY}
	BUILD=python-build
    elif [ ${PYVER} = "3.12" ]; then
	conda install -y python-build wheel "cython<3"
	pip install --no-build-isolation ${NUMPY}
	BUILD=python-build
    else
	conda install -y python-build wheel ${NUMPY}
	BUILD=python-build
    fi
    rm -rf build
    PYTHONNOUSERSITE=1 PYTHONPATH= SPACEPY_RELEASE=1 ${BUILD} -w -n -x .
    conda deactivate
    ~/miniconda/bin/conda env remove -y --name ${ENVNAME}
done
#rm -rf ~/miniconda
#rm ./Miniconda3-latest-Linux-x86_64.sh

