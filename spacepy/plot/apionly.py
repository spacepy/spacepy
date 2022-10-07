"""
Formerly used to import spacepy.plot without changing styles.

.. versionchanged:: 0.5.0
   Deprecated, as styles are no longer applied on import.
"""
import warnings

warnings.warn('Plot styles no longer applied on import, apionly deprecated'
              ' in SpacePy 0.5.0', DeprecationWarning)
