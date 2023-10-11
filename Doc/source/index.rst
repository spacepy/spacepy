.. SpacePy documentation master file, created by
   sphinx-quickstart on Tue May 31 15:38:18 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SpacePy |version| documentation
===============================

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
    @article{niehof2022spacepy,
    title={The SpacePy space science package at 12 years},
    author={Niehof, Jonathan T and Morley, Steven K and Welling, Daniel T and Larsen, Brian A},
    journal={Frontiers in Astronomy and Space Sciences},
    volume={9},
    year={2022},
    doi={10.3389/fspas.2022.1023612},
    publisher={Frontiers}
    }

To cite the code itself:
    @software{spacepy_code,
    author       = {Morley, Steven K. and Niehof, Jonathan T. and
    Welling, Daniel T. and Larsen, Brian A. and
    Brunet, Antoine and Engel, Miles A. and
    Gieseler, Jan and Haiducek, John and
    Henderson, Michael and Hendry, Aaron and
    Hirsch, Michael and Killick, Peter and
    Koller, Josef and Merrill, Asher and
    Rastatter, Lutz and Reimer, Ashton and
    Shih, Albert Y. and Stricklan, Amanda},
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


.. currentmodule:: spacepy

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
    tests
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
    ctrans
    datamanager
    datamodel
    data_assimilation
    empiricals
    igrf
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
