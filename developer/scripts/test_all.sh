#!/usr/bin/env bash

# Test lots of python/dependency versions.
# Basically superfluous to the CI, but nice to be able to test easily without
# having to put in the PR....

# ffnet install will fail if LD_LIBRARY_PATH is set!

# This might be necessary for some versions of scipy that don't have binaries
sudo aptitude install libblas-dev liblapack-dev
# download/install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda
# download/install CDF
wget https://spdf.sci.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_0/linux/cdf38_0-dist-cdf.tar.gz; tar xzf cdf38_0-dist-cdf.tar.gz; pushd cdf38_0-dist; make OS=linux ENV=gnu all; make INSTALLDIR=$HOME/cdf install; popd
source ${HOME}/cdf/bin/definitions.B

# For the most part, these versions are: what's the earliest we support,
# coupled with what's the earliest and latest version of each dep
# which works with that Python.
TESTS=(
       "2.7|numpy>=1.10.0,<1.11.0|scipy>=0.11.0,<0.12.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.0,<1.1 h5py>=2.6,<2.7 ffnet>=0.7.0,<0.8 astropy>=1.0,<1.1"
       "2.7|numpy>=1.16.0,<1.17.0|scipy matplotlib networkx h5py ffnet astropy"
       "3.5|numpy>=1.10.0,<1.11.0|scipy>=0.17.0,<0.18.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.3,<1.4 h5py>=2.6,<2.7 ffnet>=0.8.0,<0.9 astropy>=1.0,<1.1"
       "3.5|numpy>=1.18.0,<1.19.0|scipy matplotlib networkx h5py ffnet astropy"
       "3.6|numpy>=1.12.0,<1.13.0|scipy>=0.19.0,<0.20.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.3,<1.4 h5py>=2.6,<2.7 ffnet>=0.8.0,<0.9 astropy>=1.0,<1.1"
       "3.6|numpy>=1.19.0|scipy matplotlib networkx h5py ffnet astropy"
       "3.7|numpy>=1.15.1,<1.16.0|scipy>=1.0.0,<1.1.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.3,<1.4 h5py>=2.6,<2.7 ffnet>=0.8.0,<0.9 astropy>=2.0,<2.1"
       "3.7|numpy>=1.19.0|scipy matplotlib networkx h5py ffnet astropy"
       "3.8|numpy>=1.17.0,<1.18.0|scipy>=1.0.0,<1.1.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.3,<1.4 h5py>=2.6,<2.7 ffnet>=0.8.0,<0.9 astropy>=2.0,<2.1"
       "3.8|numpy>=1.19.0|scipy matplotlib networkx h5py ffnet astropy"
# numpy 1.17 works on 3.9, but doesn't build ffnet, whereas 1.18 does.
# So assuming this is an f2py problem, go for 1.18 minimum.
       "3.9|numpy>=1.18.0,<1.19.0|scipy>=1.5.0,<1.6.0 matplotlib>=1.5.0,<1.6.0 networkx>=1.3,<1.4 h5py>=2.6,<2.7 ffnet>=0.8.0,<0.9 astropy>=2.0,<2.1"
       "3.9|numpy>=1.19.0|scipy matplotlib networkx h5py ffnet astropy"
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
