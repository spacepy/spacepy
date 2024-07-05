#!/usr/bin/env bash

# This must run INSIDE docker, from the root of spacepy checkout:
# docker run -i -t --rm -v `pwd`:/io quay.io/pypa/manylinux2014_x86_64 /bin/bash /io/developer/scripts/manylinux_wheels.sh
cd /tmp
cp -a /io spacepy
/usr/bin/curl https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_1/linux/cdf39_1-dist-cdf.tar.gz --output cdf39_1-dist-cdf.tar.gz
tar -zxvf cdf39_1-dist-cdf.tar.gz
cd cdf39_1-dist
make OS=linux ENV=gnu SHARED=yes CURSES=no FORTRAN=no all
cp -a src/lib/libcdf.so /tmp/spacepy/libcdf.so
cd /tmp/spacepy

VERSIONS=(
    "cp36-cp36m"
    "cp37-cp37m"
    "cp38-cp38"
    "cp39-cp39"
    "cp310-cp310"
    "cp311-cp311"
    "cp312-cp312"
)
for PYVER in "${VERSIONS[@]}"
do
    OLDPATH=${PATH}
    PATH=/opt/python/${PYVER}/bin:${PATH}
    if [ ${PYVER} = "cp311-cp311" ]; then
	# Arbitrary version, I just know 69 doesn't work and this does...
	pip install build "setuptools<64"
    else
	pip install build setuptools
    fi
    rm -rf build
    SPACEPY_RELEASE=1 pyproject-build -w -n -x .
    pip install auditwheel
    auditwheel repair --plat manylinux2014_x86_64 dist/spacepy-*${PYVER}-linux_x86_64.whl
    PATH=${OLDPATH}
done
cp -a wheelhouse/spacepy-*manylinux_2_17_x86_64.manylinux2014_x86_64.whl /io/dist/
