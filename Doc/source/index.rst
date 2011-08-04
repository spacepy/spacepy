.. SpacePy documentation master file, created by
   sphinx-quickstart on Tue May 31 15:38:18 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SpacePy documentation
===================================

SpacePy is a package for Python, targeted at the space sciences, that aims to
make basic data analysis, modeling and visualization easier. It builds on the
capabilities of the well-known NumPy and MatPlotLib packages. Publication
quality output direct from analyses is emphasized among other goals:

    + Quickly obtain data
    + Create publications quality plots
    + Perform complicated analysis easily
    + Run common empirical models
    + Change coordinates effortlessly
    + Harness the power of Python

The SpacePy project seeks to promote accurate and open research standards by
providing an open environment for code development. In the space physics
community there has long been a significant reliance on proprietary languages
that restrict free transfer of data and reproducibility of results. By
providing a comprehensive, open-source library of widely-used analysis and
visualization tools in a free, modern and intuitive language, we hope that
this reliance will be diminished.


When publishing research which used SpacePy, please provide appropriate
credit to the SpacePy team via citation or acknowledgement.

To cite SpacePy in publications, use (BibTeX code):
@INPROCEEDINGS{spacepy11,
author = {{Morley}, S.~K. and {Koller}, J. and {Welling}, D.~T. and {Larsen}, B.~A. and {Henderson}, M.~G. and {Niehof}, J.~T.},
title = "{Spacepy - A Python-based library of tools for the space sciences}",
booktitle = "{Proceedings of the 9th Python in science conference (SciPy 2010)}",
year = 2011,
address = {Austin, TX}
}

Certain modules may provide additional citations in the __citation__
attribute. Contact a module's author before publication or public
presentation of analysis performed by that module. This allows the author
to validate the analysis and receive appropriate credit for his or her
work.


.. module:: spacepy

SpacePy Documents:

.. toctree::
    :maxdepth: 1

    doc_standard
    tips

SpacePy Code:

.. toctree::
    :maxdepth: 1

    borg
    coordinates
    datamodel
    empiricals
    omni
    poppy
    radbelt
    seapy
    time
    toolbox



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


--------------------------

:Release: |version|
:Doc generation date: |today|
