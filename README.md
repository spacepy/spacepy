[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3252523.svg)](https://doi.org/10.5281/zenodo.3252523)
[![Build Status](https://github.com/spacepy/spacepy/workflows/CI/badge.svg?branch=master)](https://github.com/spacepy/spacepy/actions?query=workflow%3ACI)

# SpacePy

SpacePy is a package for Python, targeted at the space sciences, that aims to make basic data analysis, modeling and visualization easier. It builds on the capabilities of the well-known NumPy and MatPlotLib packages. Publication quality output direct from analyses is emphasized among other goals:

 - Quickly obtain data
 - Read (and write) data from (and to) data formats like NASA CDF and HDF5
 - Create publications quality plots
 - Perform complicated analysis easily
 - Run common empirical models
 - Change coordinates and time systems effortlessly
 - Harness the power of Python

The SpacePy project seeks to promote accurate and open research standards by providing an open environment for code development. In the space physics community there has long been a significant reliance on proprietary languages that restrict free transfer of data and reproducibility of results. By providing a comprehensive, open-source library of widely-used analysis and visualization tools in a free, modern and intuitive language, we hope that this reliance will be diminished.

To help foster an open and welcoming environment, we have adopted a [code of conduct](https://github.com/spacepy/spacepy/blob/master/code-of-conduct.md) that we encourage members of the SpacePy community to read and follow.

## Getting SpacePy

Our latest release version is available through PyPI and can be installed using

```
pip install spacepy --user
```

This will also automatically install most dependencies. To permit binary installations without a compiler, this will not install ffnet on Windows. Users needing the LANLstar module can install ffnet separately (requires Fortran compiler); this can be done before or after the SpacePy install.

The latest "bleeding-edge" source code is available from our github repository at [https://github.com/spacepy/spacepy](https://github.com/spacepy/spacepy) and can be installed using the standard

```
python setup.py install --user
```

Further installation documentation can be found [here](https://spacepy.github.io/install.html) Mac-specific information can be found [here](https://spacepy.github.io/install_mac.html)
Full documentation is at [https://spacepy.github.io](https://spacepy.github.io)

SpacePy supports both Python 2.7 and 3.x.

### Dependencies

SpacePy has a number of well-maintained dependencies, most of which are automatically installed by ```pip```. These include:
 - numpy (>=1.10, !=1.15.0)
 - scipy (>=0.11)
 - matplotlib (>=1.5)
 - h5py

Soft dependencies (that are required only for a very limited part of SpacePy's functionality) are:
 - ffnet
 - NASA CDF

For complete installation, excepting pre-built Windows binaries, SpacePy also requires C and Fortran compilers. We test with GCC compilers but try to maintain support for all major compilers.

#### NASA CDF
If you wish to use CDF files, download and install the NASA CDF library. The default installation directory is recommended to help SpacePy find the library. Get the package from [https://cdf.gsfc.nasa.gov/html/sw_and_docs.html](https://cdf.gsfc.nasa.gov/html/sw_and_docs.html)

## Attribution

When publishing research which used SpacePy, please provide appropriate credit to the SpacePy team via citation or acknowledgement.

To cite SpacePy in publications, use (BibTeX code):

```
@INPROCEEDINGS{spacepy11,
author = {{Morley}, S.~K. and {Koller}, J. and {Welling}, D.~T. and {Larsen}, B.~A. and {Henderson}, M.~G. and {Niehof}, J.~T.},
title = "{Spacepy - A Python-based library of tools for the space sciences}",
booktitle = "{Proceedings of the 9th Python in science conference (SciPy 2010)}",
year = 2011,
address = {Austin, TX}
}
```

Certain modules may provide additional citations in the ```__citation__``` attribute. Contact a module's author before publication or public presentation of analysis performed by that module. This allows the author to validate the analysis and receive appropriate credit for his or her work.

For acknowledging SpacePy, please provide the URL to our github repository. [github.com/spacepy/spacepy](https://github.com/spacepy/spacepy)


