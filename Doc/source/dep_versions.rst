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
    dependencies. E.g. SciPy 1.9 requires NumPy 1.18. Although SpacePy
    supports SciPy 1.9 and Numpy 1.17, it contains no workarounds for
    using them in that combination.
 #. Support for a particular version of a dependency does not imply
    a commitment to work around bugs in that version.
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

.. list-table:: SpacePy dependency versions (2022/9/20)
   :widths: 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dependency
     - Current Release
     - Ubuntu 18.04LTS (reference only)
     - Ubuntu 20.04LTS
     - Ubuntu 22.04LTS
     - NEP 29 (42/24 mo.)
     - NEP 29 (2/3 minor versions)
     - SpacePy current minimum
   * - `CPython <https://www.python.org/downloads/>`_
     - 3.10.7 (2022/9/6)
     - 3.6.5 (2018/2/5)
     - 3.8.2 (2020/2/24)
     - 3.10.4 (2022/3/24)
     - **3.8.0** (2019/10/14)
     - 3.9.0 (2020/10/5)
     - 3.6.0 (2016/12/23)
   * - `AstroPy <https://docs.astropy.org/en/stable/changelog.html#changelog>`_
     - 5.1 (2022/5/24)
     - 2.0.4 (2018/2/6)
     - **4.0** (2019/12/16)
     - 5.0.2 (2022/3/10)
     - 4.1 (2020/10/21)
     - 4.3 (2021/7/26)
     - 1.0 (2015/2/18)
   * - `CDF <https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest-release/unix/CHANGES.txt>`_
     - 3.8.1 (2021/10/27)
     - N/A
     - N/A
     - N/A
     - 3.8.0 (2019/10/27)
     - **3.6.0** (2015/2/5)
     - 3.5.0 (2013/2/25)
   * - `dateutil <https://github.com/dateutil/dateutil/releases>`_
     - 2.8.2 (2021/7/8)
     - 2.6.1 (2017/7/10)
     - 2.7.3 (2018/5/9)
     - 2.8.1 (2019/11/3)
     - 2.8.0 (2019/2/5)
     - **2.6.0** (2016/11/8)
     - 2.1.0 (2021/3/28)
   * - `h5py <https://github.com/h5py/h5py/releases>`_
     - 3.7.0 (2022/5/24)
     - 2.7.1 (2017/9/1)
     - **2.10.0** (2019/9/6)
     - 3.6.0 (2021/11/16)
     - 3.0.0 (2020/10/30)
     - 3.5.0 (2021/10/20)
     - tested from 2.6 (2017/3/18)
   * - `matplotlib <https://github.com/matplotlib/matplotlib/releases>`_
     - 3.6.0 (2022/9/16)
     - 2.1.1 (2017/12/9)
     - **3.1.2** (2019/12/4)
     - 3.5.1 (2021/12/11)
     - 3.4.0 (2021/3/26)
     - 3.4.0 (2021/3/26)
     - 1.5.0 (2015/10/29)
   * - `numpy <https://github.com/numpy/numpy/releases>`_
     - 1.23.3 (2022/9/9)
     - 1.13.3 (2017/7/6)
     - **1.16.5** (2019/8/27)
     - 1.21.5 (2021/12/19)
     - 1.20.0 (2021/1/30)
     - 1.21.0 (2021/6/22)
     - 1.12.0 (2017/1/15)
   * - `scipy <https://github.com/scipy/scipy/releases>`_
     - 1.9.1 (2022/8/26)
     - 0.19.1 (2017/6/23)
     - **1.3.3** (2019/11/23)
     - 1.8.0 (2022/2/5)
     - 1.6.0 (2020/12/31)
     - 1.7.0 (2021/6/20)
     - 0.19.0 (2017/3/9)
   * - `sphinx <https://www.sphinx-doc.org/en/master/changes.html>`_
       (only needed for developers to build documentation)
     - 5.1.1 (2022/7/26)
     - 1.6.7 (2018/2/4)
     - **1.8.5** (2019/3/10)
     - 4.3.2 (2021/12/19)
     - 3.3.0 (2020/11/2)
     - 4.5.0 (2022/3/28)
     - tested from 1.3 (2015/3/10)
