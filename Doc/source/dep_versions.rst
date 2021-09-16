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

Regardless of minimum requirements, using the latest stable version of
a package is generally preferred. The minimum supported version for
SpacePy may not be recommended for other reasons (e.g. bug fixes or
improved features elsewhere in the package.)

This table summarizes the versions to be supported according to the
above policy, as well as the minimum version currently supported by
SpacePy. Where available, the dependency name links to its version
history. The oldest version supported according to this policy is in
**bold**.

.. list-table:: SpacePy dependency versions (2021/9/16)
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dependency
     - Current Release
     - Ubuntu 18.04LTS
     - Ubuntu 20.04LTS
     - NEP 29 (42/24 mo.)
     - NEP 29 (2/3 minor versions)
     - SpacePy current minimum
   * - `CPython <https://www.python.org/downloads/>`_
     - 3.9.7 (2021/8/30)
     - **3.6.5** (2018/2/5)
     - 3.8.2 (2020/2/24) 
     - 3.7.0 (2018/6/27)
     - 3.8.0 (2019/10/14)
     - 3.2.0 (2011/2/20) or 2.7.0 (2010/7/3)
   * - `AstroPy <https://docs.astropy.org/en/stable/changelog.html#changelog>`_
     - 4.3.1 (2021/8/11)
     - **2.0.4** (2018/2/6)
     - 4.0 (2019/2/16)
     - 3.2 (2019/6/14)
     - 4.1 (2020/10/21)
     - 1.0 (2015/2/18)
   * - `CDF <https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest-release/unix/CHANGES.txt>`_
     - 3.8.0.1 (2020/7/7)
     - N/A
     - N/A
     - 3.8.0 (2019/10/27)
     - **3.6.0** (2015/2/5)
     - 2.7.0 (1999/9/27)
   * - `dateutil <https://github.com/dateutil/dateutil/releases>`_
     - 2.8.2 (2021/7/8)
     - 2.6.1 (2017/7/10)
     - 2.7.3 (2018/5/9)
     - 2.8.0 (2019/2/5)
     - **2.6.0** (2016/11/8)
     - tested from 1.4 (2008/2/27)
   * - `ffnet <https://github.com/mrkwjc/ffnet/releases>`_
     - 0.8.4 (2018/10/28)
     - N/A
     - N/A
     - 0.8.4 (2018/10/28)
     - **0.6.0** (2007/3/22)
     - tested from 0.7 (2011/8/8)
   * - `h5py <https://github.com/h5py/h5py/releases>`_
     - 3.4.0 (2021/8/3)
     - **2.7.1** (2017/9/1)
     - 2.10.0 (2019/9/6)
     - 3.0.0 (2020/10/30)
     - 3.2.0 (2021/3/3)
     - tested from 2.6 (2017/3/18)
   * - `matplotlib <https://github.com/matplotlib/matplotlib/releases>`_
     - 3.4.3 (2021/8/12)
     - **2.1.1** (2017/12/9)
     - 3.1.2 (2019/12/4)
     - 3.2.0 (2020/3/3)
     - 3.2.0 (2020/3/3)
     - 1.5.0 (2015/10/29)
   * - `networkx <https://github.com/networkx/networkx/releases>`_
     - 2.6.3 (2021/9/9)
     - **1.11** (2016/1/30)
     - 2.4 (2019/10/16)
     - 2.4 (2019/10/16)
     - 2.4 (2019/10/16)
     - tested from 1.0 (2010/1/7)
   * - `numpy <https://github.com/numpy/numpy/releases>`_
     - 1.21.2 (2021/8/15)
     - **1.13.3** (2017/7/6)
     - 1.16.5 (2019/8/27)
     - 1.18.0 (2019/12/22)
     - 1.19.0 (2020/1/20)
     - 1.10.0 (2015/10/5)
   * - `scipy <https://github.com/scipy/scipy/releases>`_
     - 1.7.1 (2021/8/1)
     - **0.19.1** (2017/6/23)
     - 1.3.3 (2019/11/23)
     - 1.4.0 (2019/12/16)
     - 1.5.0 (2020/6/21)
     - 0.11.0 (2012/9/24)
