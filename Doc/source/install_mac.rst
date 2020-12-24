******************
MacOS Installation
******************

The following are performed from the command line on OSX.

Installation on OSX requires a C compiler::

   xcode-select --install

Our recommended (but not required) Python distribution is `Anaconda
<https://docs.anaconda.com/anaconda/>`_ running 64-bit
Python 3. Anaconda includes much of the scientific Python
stack. Another excellent distribution is `Canopy
<https://www.enthought.com/product/canopy/>`_.

Once python is installed (Anaconda assumed) install the Fortran compiler::

   conda install gfortran_osx-64


Dependencies via conda and pip
==============================

Installation via ``pip`` will automatically install most Python
dependencies (but not the :ref:`NASA CDF library <linux_CDF>`). They
can also be installed from conda::

   conda install numpy scipy matplotlib networkx h5py
   pip install ffnet

Once this is set up, ``pip install spacepy`` should just work. If
you're installing as a single user (not in a virtual environment) then
add the ``--user`` flag.

You will also need the :ref:`NASA CDF library <linux_CDF>` to use
:mod:`~spacepy.pycdf`.

``ffnet`` is required for :mod:`~spacepy.LANLstar`. It is
automatically installed when SpacePy is installed from source
(currently the only supported means on Mac). It can also be explicitly
installed via pip::

  pip install ffnet

To install the latest code from the repository, rather than
the latest stable release, use::

   pip install git+https://github.com/spacepy/spacepy

.. contents::
   :local:


