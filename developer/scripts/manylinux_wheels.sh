#!/usr/bin/env bash

# This must run INSIDE docker, from the root of spacepy checkout:
# docker run -i -t --rm -v `pwd`:/io quay.io/pypa/manylinux2014_x86_64 /bin/bash /io/developer/scripts/manylinux_wheels.sh
# for ARM:
# docker run -i -t --rm -v `pwd`:/io --platform linux/arm64 quay.io/pypa/manylinux2014_aarch64 /bin/bash /io/developer/scripts/manylinux_wheels.sh
cd /tmp
cp -a /io spacepy
/usr/bin/curl https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_1/linux/cdf39_1-dist-cdf.tar.gz --output cdf39_1-dist-cdf.tar.gz
tar -zxvf cdf39_1-dist-cdf.tar.gz
cd cdf39_1-dist
make OS=linux ENV=gnu SHARED=yes CURSES=no FORTRAN=no all
cp -a src/lib/libcdf.so /tmp/spacepy/libcdf.so
cd /tmp/spacepy

ARCH=`uname -m`
VERSIONS=(
    "cp310-cp310"
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
    auditwheel repair --plat manylinux2014_${ARCH} dist/spacepy-*-cp36-abi3-linux_${ARCH}.whl
    PATH=${OLDPATH}
done
mkdir -p /io/dist/
cp -a wheelhouse/spacepy-*manylinux_2_17_${ARCH}.manylinux2014_${ARCH}.whl /io/dist/
