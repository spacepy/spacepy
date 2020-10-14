

############################################################
toolbox - Toolbox of various functions and generic utilities
############################################################

.. automodule:: spacepy.toolbox

.. currentmodule:: spacepy.toolbox

- `Array binning`_
- `Array creation`_
- `Array searching and masking`_
- `Other functions`_
- `Multithreading and multiprocessing`_
- `System tools`_
- `Quaternion math`_

Array binning
-------------
.. autosummary::
    :toctree: autosummary

    arraybin
    bin_center_to_edges
    bin_edges_to_center
    binHisto

Array creation
--------------
.. autosummary::
    :toctree: autosummary

    dist_to_list
    geomspace
    linspace
    logspace

Array searching and masking
---------------------------
.. autosummary::
    :toctree: autosummary

    interweave
    isview
    tCommon
    tOverlap
    tOverlapHalf

Other functions
---------------
.. autosummary::
    :toctree: autosummary

    assemble
    bootHisto
    dictree
    eventTimer
    feq
    getNamedPath
    human_sort
    hypot
    indsFromXrange
    interpol
    intsolve
    medAbsDev
    mlt2rad
    normalize
    pmm
    poisson_fit
    rad2mlt
    windowMean

Multithreading and multiprocessing
----------------------------------
.. autosummary::
    :toctree: autosummary

    thread_job
    thread_map

System tools
------------
.. autosummary::
    :toctree: autosummary

    do_with_timeout
    get_url
    loadpickle
    progressbar
    query_yes_no
    savepickle
    timeout_check_call
    TimeoutError
    update

.. _toolbox_quaternions:

Quaternion math
---------------

Quaternion functions have moved to :mod:`~spacepy.coordinates`;
new functions are also added there.

.. autosummary::
    :toctree: autosummary

    quaternionRotateVector
    quaternionNormalize
    quaternionMultiply
    quaternionConjugate
