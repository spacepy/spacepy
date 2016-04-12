"""
Import spacepy.plot without changing styles.

This will import :mod:`spacepy.plot` without applying the default SpacePy
plot styles. Access all functions and classes through attributes of
 :mod:`~spacepy.plot`; nothing further is defined here.
"""

from . import revert_style
revert_style()
del revert_style
