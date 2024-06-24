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
       *not* necessarily be tested in continuous integration. In particular,
       CI will only be run against versions of dependencies other than numpy
       which have binary wheels available. Numpy may be tested with earlier
       versions that do not conflict with other dependencies. For this reason,
       the CI configuration is not a reasonable guide of the minimum supported
       versoin.

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

.. list-table:: SpacePy dependency versions (2024/6/24)
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dependency
     - Current Release
     - Ubuntu 24.04LTS
     - Ubuntu 22.04LTS
     - NEP 29 (42/24 mo.)
     - NEP 29 (2/3 minor versions)
     - SpacePy current minimum
   * - `CPython <https://www.python.org/downloads/>`_
     - 3.12.4 (2024/6/6)
     - 3.12.3 (2024/4/9)
     - 3.10.4 (2022/3/24)
     - **3.10.0** (2021/10/4)
     - 3.11.0 (2022/10/24)
     - 3.6.0 (2016/12/23)
   * - `AstroPy <https://docs.astropy.org/en/stable/changelog.html#changelog>`_
     - 6.1.1 (2024/6/14)
     - 6.0.0 (2023/11/25)
     - **5.0.2** (2022/3/10)
     - 5.2 (2022/12/12)
     - 5.3 (2023/5/22)
     - 1.0 (2015/2/18)
   * - `CDF <https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest-release/unix/CHANGES.txt>`_
     - 3.9.0 (2023/1/22)
     - N/A
     - N/A
     - 3.9.0 (2023/1/22)
     - **3.7.0** (2018/5/11)
     - 3.5.0 (2013/2/25)
   * - `dateutil <https://github.com/dateutil/dateutil/releases>`_
     - 2.9.0 (2024/2/29)
     - 2.8.2 (2021/7/8)
     - 2.8.1 (2019/11/3)
     - 2.9.0 (2024/2/29)
     - **2.7.0** (2018/3/11)
     - 2.5.0 (2016/2/28)
   * - `h5py <https://github.com/h5py/h5py/releases>`_
     - 3.11.0 (2024/4/10)
     - 3.10.0 (2023/10/9)
     - **3.6.0** (2021/11/16)
     - 3.8.0 (2023/1/23)
     - 3.9.0 (2023/6/20)
     - 2.10.0 (2019/9/6)
   * - `matplotlib <https://github.com/matplotlib/matplotlib/releases>`_
     - 3.9.0 (2024/5/15)
     - 3.6.3 (2023/1/11)
     - **3.5.1** (2021/12/11)
     - 3.6.0 (2022/9/16)
     - 3.7.0 (2023/2/13)
     - 3.1.0 (2019/5/18)
   * - `numpy <https://github.com/numpy/numpy/releases>`_
     - 2.0.0 (2024/6/16)
     - 1.26.4 (2024/2/5)
     - **1.21.5** (2021/12/19)
     - 1.24.0 (2022/12/18)
     - 1.25.0 (2023/6/17)
     - 1.15.1 (2018/8/21)
   * - `pandas <https://pandas.pydata.org/docs/whatsnew/>`_
     - 2.2.2 (2024/4/10)
     - 2.1.4 (2023/12/8)
     - **1.3.5** (2022/3/10)
     - 1.5.0 (2022/9/19)
     - 2.0.0 (2023/4/3)
     - 0.18.0 (2016/3/13)
   * - `scipy <https://github.com/scipy/scipy/releases>`_
     - 1.13.1 (2024/5/23)
     - 1.11.4 (2023/11/18)
     - **1.8.0** (2022/2/5)
     - 1.9.0 (2022/7/29)
     - 1.11.0 (2023/6/25)
     - 1.0.0 (2017/10/25)
   * - `pip <https://pip.pypa.io/en/stable/news/>`_
     - 24.1 (2024/6/20)
     - 24.0 (2024/2/3)
     - **22.0.2** (2022/1/30)
     - 22.2 (2022/7/21)
     - 23.3 (2023/10/15)
     - tested with 20.0.2
   * - `setuptools <https://setuptools.pypa.io/en/latest/history.html>`_
     - 70.1.0 (2024/6/19)
     - 68.1.2 (2023/8/18)
     - **59.6.0** (2021/12/12)
     - 63.0.0 (2022/7/4)
     - 69.4.0 (2024/4/12)
     - tested with 44.1.1
   * - `wheel <https://wheel.readthedocs.io/en/stable/news.html>`_
     - 0.43.0 (2024/3/11)
     - 0.42.0 (2023/11/26)
     - **0.37.1** (2021/12/22)
     - 0.38.0 (2022/10/21)
     - 0.41.0 (2023/8/5)
     - tested with 0.34.2
   * - `sphinx <https://www.sphinx-doc.org/en/master/changes.html>`_
       (only needed for developers to build documentation)
     - 7.3.7 (2024/4/19)
     - 7.2.6 (2023/9/13)
     - **4.3.2** (2021/12/19)
     - 5.1.0 (2022/7/24)
     - 7.1.0 (2023/7/24)
     - 1.8.0 (2018/9/13)
   * - `build <https://pypa-build.readthedocs.io/en/latest/changelog.html>`_
       (only needed for developers to build releases)
     - 1.2.1 (2024/3/28)
     - 1.0.3 (2023/9/6)
     - **0.7.0** (2021/9/16)
     - 0.9.0 (2022/10/27)
     - 1.0.0 (2023/9/1)
     - tested with 0.4.0
