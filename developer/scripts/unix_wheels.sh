#!/usr/bin/env bash

# Get miniconda if you need it
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda
VERSIONS=(
    "2.7|numpy>=1.10.0,<1.11.0"
    "3.5|numpy>=1.10.0,<1.11.0"
    "3.6|numpy>=1.12.0,<1.13.0"
    "3.7|numpy>=1.15.1,<1.16.0"
    "3.8|numpy>=1.17.0,<1.18.0"
# numpy 1.17 works on 3.9, but doesn't build ffnet, whereas 1.18 does.
# So assuming this is an f2py problem, go for 1.18 minimum.
    "3.9|numpy>=1.18.0,<1.19.0"
)
for thisBuild in "${VERSIONS[@]}"
do
    IFS='|' read -ra tmpVar <<< "$thisBuild"
    PYVER=${tmpVar[0]}
    NUMPY="${tmpVar[1]}"
    ENVNAME=spacepy${PYVER//.}
    ~/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
    source ~/miniconda/bin/activate ${ENVNAME}
    pip install ${NUMPY}
    rm -rf build
    PYTHONNOUSERSITE=1 PYTHONPATH= python setup.py bdist_wheel
    conda deactivate
    ~/miniconda/bin/conda env remove --name ${ENVNAME}
done
#rm -rf ~/miniconda
#rm ./Miniconda3-latest-Linux-x86_64.sh

