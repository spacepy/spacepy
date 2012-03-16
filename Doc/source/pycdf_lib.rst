#################################
pycdf.lib - Access to CDF library
#################################

.. currentmodule:: spacepy.pycdf
.. data:: lib

Module global :class:`Library` object.

Initalized at :mod:`~spacepy.pycdf` load time so all classes have ready
access to the CDF library and a common state. E.g:

>>> from spacepy import pycdf
>>> pycdf.lib.version
    (3, 3, 0, ' ')