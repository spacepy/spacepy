#!/usr/bin/env bash
#bash script to make SpacePy tarball

/usr/bin/env python setup.py sdist --formats=gztar,zip
