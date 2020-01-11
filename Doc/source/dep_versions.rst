**************************
Dependency version support
**************************

SpacePy will occasionally drop support for old versions of
:doc:`dependencies <dependencies>`. Failures with older versions will
not be treated as SpacePy bugs. Dependency support is based on these
principles:

 #. SpacePy supports released versions of a dependency that meet a
    minimum version requirement; there is no maximum supported
    version.
 #. Support for old versions of dependencies will be dropped only for
    reason, e.g. if a new version is required to support a new feature
    or fix a bug. Maintenance of convoluted workarounds is included in
    this category.
 #. Support will be maintained for versions included in the
    second-most-recent Ubuntu Long Term Support (LTS) release,
    e.g. upon release of Ubuntu 20.04 LTS, support will be maintained
    for at least the versions in 18.04 LTS
 #. Support will also be maintained for Python and NumPy versions
    in the spirit of `NEP 29
    <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_.

    #. A SpacePy release will support at least all minor versions of Python
       released in the prior 42 months, and at least the two latest minor
       versions.
    #. SpacePy will support all minor versions of NumPy and, where
       possible, other Python dependencies released in the prior 24 months,
       and at least the three latest minor versions.
    #. Non-Python dependencies that use a similar versioning system will
       be supported similarly where possible.
    #. This support is based on minor releases (x.y.0), not subsequent
       subminor releases (x.y.z for the same x.y). Where x.y.0 is supported,
       so is x.y.z for all z.
    #. All versions of all dependencies and all combinations thereof will
       *not* necessarily be tested in continuous integration.

 #. No support will be provided for conflicting versions of
    dependencies. E.g. SciPy 1.4 requires NumPy 1.13. Although SpacePy
    supports SciPy 1.4 and Numpy 1.6, it contains no workarounds for
    using them in that combination.
 #. Support for a particular version of a dependency does not imply
    a commitment to work around bugs in that version.
 #. Support for even earlier versions will be maintained as necessary
    for Python 2 compatibility as long as :doc:`Python 2 support
    <py2k_eol>` is maintained.
 #. A release of SpacePy that requires new dependency versions will
    always have a subminor version of 0, e.g. if the release that
    follows 0.5.2 requires updated dependencies, it will be numbered
    0.6.0.
 #. The commit that requires a newer version of a dependency must also
    update the ``requirements.txt``, :doc:`dependencies`, and the
    table below. The commit message must include the reason for the
    dependency requirement.

This table summarizes the versions to be supported according to the
above policy, as well as the minimum version currently supported by
SpacePy. Where available, the dependency name links to its version
history. The oldest version supported according to this policy is in
**bold**.

.. list-table:: SpacePy dependency versions (2020/1/9)
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dependency
     - Current Release
     - Ubuntu 16.04LTS
     - Ubuntu 18.04LTS
     - NEP 29 (42/24 mo.)
     - NEP 29 (2/3 minor versions)
     - SpacePy current minimum
   * - `CPython <https://www.python.org/downloads/>`_
     - 3.8.1 (2019/12/18)
     - **3.5.1** (2015/12/7)
     - 3.6.5 (2018/2/5)
     - 3.6.0 (2016/12/23)
     - 3.7.0 (2018/6/27)
     - 3.2.0 (2011/2/20) or 2.7.0 (2010/7/3)
   * - `CDF <https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest-release/unix/CHANGES.txt>`_
     - 3.7.1 (2018/8/21)
     - N/A
     - N/A
     - 3.7.0 (2018/5/11)
     - **3.5.0** (2013/9/15)
     - 2.7.0 (1999/9/27)
   * - `dateutil <https://github.com/dateutil/dateutil/releases>`_
     - 2.8.1 (2019/11/3)
     - **2.4.2** (2015/3/31)
     - 2.6.1 (2017/7/10)
     - 2.7.0 (2018/3/11)
     - 2.6.0 (2016/11/8)
     - minimum not established
   * - `ffnet <https://github.com/mrkwjc/ffnet/releases>`_
     - 0.8.4 (2018/10/28)
     - N/A
     - N/A
     - 0.8.4 (2018/10/28)
     - **0.6.0** (2007/3/22)
     - minimum not established (tested from 0.7.0)
   * - `h5py <https://github.com/h5py/h5py/releases>`_
     - 2.10.0 (2019/9/6)
     - **2.6.0** (2017/3/18)
     - 2.7.1 (2017/9/1)
     - 2.8.0 (2018/5/13)
     - 2.8.0 (2018/5/13)
     - minimum not established
   * - `matplotlib <https://github.com/matplotlib/matplotlib/releases>`_
     - 3.1.2 (2019/12/4)
     - **1.5.1** (2016/1/10)
     - 2.1.1 (2017/12/9)
     - 2.1.0 (2017/10/7)
     - 3.0.0 (2018/9/17)
     - 1.5.0 (2015/10/29)
   * - `networkx <https://github.com/networkx/networkx/releases>`_
     - 2.4.0 (2019/10/16)
     - 1.11 (2016/1/30)
     - **1.11** (2016/1/30)
     - 2.1 (2018/1/22)
     - 2.1 (2018/1/22)
     - minimum not established
   * - `numpy <https://github.com/numpy/numpy/releases>`_
     - 1.18.1 (2020/1/6)
     - **1.11.0** (2016/3/27)
     - 1.13.1 (2017/7/6)
     - 1.14.0 (2018/1/6)
     - 1.16.0 (2019/1/13)
     - 1.4.0 (2009/12/27)
   * - `scipy <https://github.com/scipy/scipy/releases>`_
     - 1.4.1 (2019/12/18)
     - **0.17.0** (2016/1/23)
     - 0.19.1 (2017/6/23)
     - 1.2.0 (2018/12/17)
     - 1.1.0 (2018/5/5)
     - 0.10.0 (2011/11/13)
