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
credit to the SpacePy team via citation or acknowledgment.

To cite SpacePy in publications, use (BibTeX code):
    @INPROCEEDINGS{spacepy11,
    author = {{Morley}, S.~K. and {Koller}, J. and {Welling}, D.~T. and {Larsen}, B.~A. and {Henderson}, M.~G. and {Niehof}, J.~T.},
    title = "{Spacepy - A Python-based library of tools for the space sciences}",
    booktitle = "{Proceedings of the 9th Python in science conference (SciPy 2010)}",
    year = 2011,
    address = {Austin, TX}
    }

Or to cite the code itself:
    @software{SpacePy,
    author       = {{Larsen}, B.~A. and {Morley}, S.~K. and {Niehof}, J.~T. and {Welling}, D.~T.},
    title        = {SpacePy},
    publisher    = {Zenodo},
    doi          = {10.5281/zenodo.3252523},
    url          = {https://doi.org/10.5281/zenodo.3252523}
    }

Certain modules may provide additional citations in the ``__citation__``
attribute. Contact a module's author (details in the ``__citation__`` attribute) 
before publication or public presentation of analysis performed by that 
module, or in case of questions about the module. This allows the author to 
validate the analysis and receive appropriate credit for his or her
work.


.. module:: spacepy

Getting Started
===============

First steps in SpacePy and scientific Python.

.. toctree::
    :maxdepth: 1

    install
    quickstart
    help

SpacePy Documents
=================

Further reference material on how to use SpacePy, and examples.

.. toctree::
    :maxdepth: 1

    capabilities
    release_notes
    case_studies
    publications
    py2k_eol
    configuration
    scripts

Developer Guide
===============

For those developing SpacePy, plus tips for all Python developers.

.. toctree::
    :maxdepth: 1

    pythonic
    tips
    dep_versions
    doc_standard
    ci

.. _module_reference:

SpacePy Module Reference
========================

Description of all functions within SpacePy, by module.

.. toctree::
    :maxdepth: 1

    spacepy
    ae9ap9
    coordinates
    datamanager
    datamodel
    data_assimilation
    empiricals
    irbempy
    lanlstar
    omni
    plot
    poppy
    pybats
    pycdf
    radbelt
    seapy
    time
    toolbox
..  realtime


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


--------------------------

:Release: |version|
:Doc generation date: |today|
