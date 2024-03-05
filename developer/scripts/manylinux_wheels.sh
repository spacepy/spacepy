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
    "cp36-cp36m|numpy~=1.12.0"
    "cp37-cp37m|numpy~=1.15.1"
    "cp38-cp38|numpy~=1.17.0"
    "cp39-cp39|numpy~=1.17.0"
    "cp310-cp310|numpy~=1.21.0"
    "cp311-cp311|numpy~=1.22.0"
    "cp312-cp312|numpy~=1.24.0"
)
for thisBuild in "${VERSIONS[@]}"
do
    IFS='|' read -ra tmpVar <<< "$thisBuild"
    PYVER=${tmpVar[0]}
    NUMPY="${tmpVar[1]}"
    OLDPATH=${PATH}
    PATH=/opt/python/${PYVER}/bin:${PATH}
    if [ ${PYVER} = "cp311-cp311" ]; then
	# Arbitrary version, I just know 69 doesn't work and this does...
	pip install build "setuptools<64" "cython<3"
    else
	pip install build setuptools "cython<3"
    fi
    pip install --no-build-isolation ${NUMPY}
    rm -rf build
    SPACEPY_RELEASE=1 pyproject-build -w -n -x .
    pip install auditwheel
    auditwheel repair --plat manylinux2014_x86_64 dist/spacepy-*${PYVER}-linux_x86_64.whl
    PATH=${OLDPATH}
done
cp -a wheelhouse/spacepy-*manylinux_2_17_x86_64.manylinux2014_x86_64.whl /io/dist/
