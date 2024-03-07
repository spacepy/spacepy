[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3252523.svg)](https://doi.org/10.5281/zenodo.3252523)
[![Build Status](https://github.com/spacepy/spacepy/workflows/CI/badge.svg?branch=main)](https://github.com/spacepy/spacepy/actions?query=workflow%3ACI)

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

To help foster an open and welcoming environment, we have adopted a [code of conduct](https://github.com/spacepy/spacepy/blob/main/code-of-conduct.md) that we encourage members of the SpacePy community to read and follow.

## Getting SpacePy

Our latest release version is available through PyPI and can be installed using

```
pip install spacepy --user
```

This will also automatically install most dependencies.

The latest "bleeding-edge" source code is available from our github repository at [https://github.com/spacepy/spacepy](https://github.com/spacepy/spacepy).

Further installation documentation, including building from source and OS-specific information, can be found [here](https://spacepy.github.io/install.html). Full documentation is at [https://spacepy.github.io](https://spacepy.github.io).

SpacePy supports Python 3.6 and later.

### Dependencies

SpacePy has a number of well-maintained dependencies which are automatically installed by ```pip```. These include:

 - numpy (>=1.15.1)
 - dateutil (>=2.1)
 - scipy (>=1.0)
 - matplotlib (>=3.1)
 - h5py (>=2.10)

## Attribution

When publishing research which used SpacePy, please provide appropriate credit to the SpacePy team via citation or acknowledgement.

To cite SpacePy in publications, please cite both the code (DOI: 10.5281/zenodo.3252523) and the papers describing the package (BibTeX code):

```
@article{niehof2022spacepy,
  title={The SpacePy space science package at 12 years},
  author={Niehof, Jonathan T and Morley, Steven K and Welling, Daniel T and Larsen, Brian A},
  journal={Frontiers in Astronomy and Space Sciences},
  volume={9},
  year={2022},
  doi={10.3389/fspas.2022.1023612},
  publisher={Frontiers}
}
```

and/or

```
@INPROCEEDINGS{spacepy11,
author = {{Morley}, S.~K. and {Koller}, J. and {Welling}, D.~T. and {Larsen}, B.~A. and {Henderson}, M.~G. and {Niehof}, J.~T.},
title = "{Spacepy - A Python-based library of tools for the space sciences}",
booktitle = "{Proceedings of the 9th Python in science conference (SciPy 2010)}",
year = 2011,
address = {Austin, TX}
}
```

For additional information, see the [CITATION.cff](https://github.com/spacepy/spacepy/blob/main/CITATION.cff) file.
Certain modules may provide additional citations in the ```__citation__``` attribute. Contact a module's author before publication or public presentation of analysis performed by that module. This allows the author to validate the analysis and receive appropriate credit for his or her work.

For acknowledging SpacePy, please provide the URL to our github repository: [github.com/spacepy/spacepy](https://github.com/spacepy/spacepy).

## Changes
Changes in the released version of SpacePy are provided in the [release notes](https://spacepy.github.io/release_notes.html). For changes since the latest release, see the [repository version](https://github.com/spacepy/spacepy/blob/main/Doc/source/release_notes.rst).
