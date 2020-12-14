=======================
Python 2 End of Support
=======================

Python 3 is fully supported by SpacePy and used in daily work by most
of the SpacePy team.

On January 1, 2020, Python 2 reached `end of life
<https://www.python.org/doc/sunset-python-2/>`_. Most of the
scientific Python stack has committed to ending Python 2 support `by
2020 <https://python3statement.org/>`_. This includes `numpy
<https://numpy.org/neps/nep-0014-dropping-python2.7-proposal.html>`_.

As a result, the SpacePy team will phase out Python 2 support over the
course of 2020 and early 2021. This process is managed through SpacePy `issue 26
<https://github.com/spacepy/spacepy/issues/26>`_.

0.2 series: full support
========================

The last release of the SpacePy 0.2 series will be 0.2.3, by the end
of 2020. This will be the last release where all
functionality works with Python 2 and that has binary installers
provided for Python 2. 0.2 will *not* be supported past this release
(this will not be a "long-term support" release.)

0.3 series: no feature support
==============================

Starting with 0.3.0 in early 2021, the SpacePy team will:

 * Provide no prebuilt packages for Python 2. We will attempt to
   ensure the last 0.2.x version will still install from pip on
   Python 2.
 * Allow new features that do not support Python 2 as long as they do
   not break existing functionality.
 * Provide no workarounds for dependencies that no longer support
   Python 2.

SpacePy 0.3.x will still function on Python 2 for those who install
"by hand".

0.4 series: no bugfix support
=============================   

Starting with 0.4.0, no later than mid-2021, the SpacePy team will
provide no fixes for bugs that cannot be reproduced on Python 3.

SpacePy 0.4.x will still function on Python 2 for those who install
"by hand".

0.5 series: remove support
==========================

Starting with 0.5.0, mid-2021, SpacePy developers will begin
removing code that exists only to support Python 2. SpacePy 0.5.x will
not function on Python 2.

--------------------------

:Release: |version|
:Doc generation date: |today|
