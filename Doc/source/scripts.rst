***************
SpacePy Scripts
***************

Some scripts using SpacePy are included in the ``scripts`` directory
of the source distribution. At the moment they are not installed by
the installer.

istp_checks.py
==============
.. program:: istp_checks.py

Checks for various ISTP compliance issues in a file and prints any
issues found. This script is supplemental to the checker included with
the `ISTP skeleton editor <https://spdf.gsfc.nasa.gov/skteditor/>`_;
it primarily checks for errors that the skeleton editor does not.

This is just a thin wrapper to :meth:`spacepy.pycdf.istp.FileChecks.all`.

.. option:: cdffile

    Name of the CDF file to check (required)


