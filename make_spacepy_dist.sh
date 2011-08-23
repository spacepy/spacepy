#!/usr/bin/env bash
#bash script to make SpacePy tarball

PYTHON='/usr/bin/env python'


#make documentation
cd Doc
make html

cd ..
PYTHON setup.py sdist --formats=gztar,zip

# what I like about the script is that we can add in rpm or whatever we want in here
# without having to remeber everything
