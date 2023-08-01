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
       "3.6|pip~=20.0.0 setuptools~=44.1.0 wheel~=0.34.2 cython|numpy~=1.15.1|python-dateutil~=2.1.0 scipy~=1.0.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=1.0,<1.1|old"
       "3.6|pip setuptools wheel cython|numpy>=1.19.0,<1.20.0|python-dateutil scipy matplotlib h5py astropy|new"
       "3.7|pip~=20.0.0 setuptools~=44.1.0 wheel~=0.34.2 cython|numpy>=1.15.1,<1.16.0|python-dateutil~=2.1.0 scipy>=1.0.0,<1.1.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1|old"
       "3.7|pip setuptools wheel cython|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy|new"
       "3.8|pip~=20.0.0 setuptools~=44.1.0 wheel~=0.34.2 cython|numpy>=1.17.0,<1.18.0|python-dateutil~=2.1.0 scipy>=1.0.0,<1.1.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1|old"
       "3.8|pip setuptools wheel cython|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy|new"
       "3.9|pip~=20.0.0 setuptools~=44.1.0 wheel~=0.34.2 cython<3.0|numpy>=1.18.0,<1.19.0|python-dateutil~=2.1.0 scipy>=1.5.0,<1.6.0 matplotlib~=3.1.0 h5py~=2.10.0 astropy>=2.0,<2.1|old"
       "3.9|pip setuptools wheel|numpy>=1.21.0 cython|python-dateutil scipy matplotlib h5py astropy|new"
       "3.10|pip~=20.0.0 setuptools~=44.1.0 wheel~=0.34.2 cython|numpy>=1.21.0,<1.22.0|python-dateutil~=2.7.0 scipy>=1.7.2,<1.8.0 matplotlib~=3.1.0 h5py>=3.6,<3.7 astropy>=4.0,<4.1|old"
       "3.10|pip setuptools wheel cython|numpy>=1.21.0|python-dateutil scipy matplotlib h5py astropy|new"
      )
for thisTest in "${TESTS[@]}"
do
    IFS='|' read -ra tmpTest <<< "$thisTest"
    PYVER=${tmpTest[0]}
    PIP="${tmpTest[1]}"
    NUMPY="${tmpTest[2]}"
    PIPLIST="${tmpTest[3]}"
    DEPSTRAT="${tmpTest[4]}"
    ENVNAME=py${PYVER//.}
    ~/miniconda/bin/conda create -y -n ${ENVNAME} python=${PYVER}
    source ~/miniconda/bin/activate ${ENVNAME}
    # Get the preferred version of install tools, numpy's dependencies
    pip install --no-build-isolation ${PIP}
    # Get numpy in first (lots uses distutils)
    pip install --no-build-isolation ${NUMPY}
    # Don't let it grab a different version of numpy later,
    # and always build against the installed numpy
    if [ "$DEPSTRAT" = "old" -a \( ${PYVER} = "3.6" -o ${PYVER} = "3.7" -o ${PYVER} = "3.8" \) ]; then
	# scipy < 1.5 doesn't build on gcc 10+ without this argument
	FOPT="-fallow-argument-mismatch" pip install --no-build-isolation ${NUMPY} ${PIPLIST}
    else
	pip install --no-build-isolation ${NUMPY} ${PIPLIST}
    fi
    rm -rf build
    # Make sure not building against user's version of libraries
    # pip install . does not cache the wheel (phew)
    PYTHONNOUSERSITE=1 PYTHONPATH= pip install --no-build-isolation --no-deps .
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
