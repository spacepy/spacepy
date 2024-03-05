#!/usr/bin/env bash

# Run unit tests on all the wheels

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
    pip install --find-links ../dist/ --only-binary spacepy spacepy
    python test_all.py > test_output_${PYVER}.txt 2>&1
    conda deactivate
    ~/miniconda/bin/conda env remove -y --name ${ENVNAME}
done
