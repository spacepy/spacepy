#!/usr/bin/env bash
#bash script to make SpacePy tarball

PYTHON='/usr/bin/env python'


#make documentation
cd Doc
make html


## TODO THIS is a hack for now
mv build/html .


cd ..
PYTHON setup.py sdist --formats=gztar,zip
