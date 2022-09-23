#!/usr/bin/env bash

# Test lots of python/dependency versions.
# Basically superfluous to the CI, but nice to be able to test easily without
# having to put in the PR....

# This might be necessary for some versions of scipy that don't have binaries
sudo aptitude install libblas-dev liblapack-dev
# download/install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda
# download/install CDF
wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/linux/cdf38_1-dist-cdf.tar.gz; tar xzf cdf38_1-dist-cdf.tar.gz; pushd cdf38_1-dist; make OS=linux ENV=gnu all; make INSTALLDIR=$HOME/cdf install; popd
source ${HOME}/cdf/bin/definitions.B

# For the most part, these versions are: what's the earliest we support,
# coupled with what's the earliest and latest version of each dep
# which works with that Python.
TESTS=(
       "3.6|numpy~=1.15.1|python-dateutil~=2.1.0 scipy~=1.0.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=1.0,<1.1"
       "3.6|numpy>=1.19.0,<1.20.0|python-dateutil scipy matplotlib h5py astropy"
       "3.7|numpy>=1.15.1,<1.16.0|python-dateutil~=2.1.0 scipy>=1.0.0,<1.1.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1"
       "3.7|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy"
       "3.8|numpy>=1.17.0,<1.18.0|python-dateutil~=2.1.0 scipy>=1.0.0,<1.1.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1"
       "3.8|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy"
       "3.9|numpy>=1.18.0,<1.19.0|python-dateutil~=2.1.0 scipy>=1.5.0,<1.6.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1"
       "3.9|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy"
       "3.10|numpy>=1.21.0,<1.22.0|python-dateutil~=2.7.0 scipy>=1.7.2,<1.8.0 matplotlib~=3.1.0  h5py>=3.6,<3.7 astropy>=4.0,<4.1"
       "3.10|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy"
      )
for thisTest in "${TESTS[@]}"
do
    IFS='|' read -ra tmpTest <<< "$thisTest"
    PYVER=${tmpTest[0]}
    NUMPY="${tmpTest[1]}"
    PIPLIST="${tmpTest[2]}"
    ENVNAME=py${PYVER//.}
    ~/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
    source ~/miniconda/bin/activate ${ENVNAME}
    # Get numpy in first (lots uses distutils)
    pip install ${NUMPY}
    # But don't let it grab a different version of numpy later
    pip install ${NUMPY} ${PIPLIST}
    rm -rf build
    # Make sure not building against user's version of libraries
    PYTHONNOUSERSITE=1 PYTHONPATH= python setup.py build
    pushd tests
    echo Python ${PYVER}
    echo ${NUMPY} ${PIPLIST}
    PYTHONNOUSERSITE=1 PYTHONPATH= python test_all.py
    popd
    conda deactivate
    ~/miniconda/bin/conda env remove --name ${ENVNAME}
done
rm -rf ~/miniconda
rm ./Miniconda3-latest-Linux-x86_64.sh
