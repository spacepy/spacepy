#!/usr/bin/env python

"""Helper functions for SpacePy unit tests

This module is importable from test scripts that are in this directory.
"""

import os.path
import re
import sys
import sysconfig
import warnings


testsdir = os.path.dirname(os.path.abspath(__file__))
"""Directory containing the unit test scripts"""

datadir = os.path.join(testsdir, 'data')
"""Directory containing unit test data"""


def add_build_to_path():
    """Adds the python build directory to the search path.

    Locates the build directory in the same repository as this test module
    and adds the (version-specific) library directories to the Python
    module search path, so the unit tests can be run against the built
    instead of installed version.

    This is run on import of this module.
    """
    # Prioritize version-specific path; py2 tends to be version-specific
    # and py3 tends to use just "lib". But only use first-matching.
    for pth in ('lib',  # Prepending, so add low-priority paths first.
                'lib.{0}-{1}.{2}'.format(sysconfig.get_platform(),
                                         *sys.version_info[:2]),
                ):
        buildpath = os.path.abspath(os.path.join(testsdir, '..', 'build', pth))
        if os.path.isdir(buildpath):
            if not buildpath in sys.path:
                sys.path.insert(0, buildpath)
            break


class assertWarns(warnings.catch_warnings):
    """Tests that a warning is raised.

    Use as a context manager.  Code within the context manager block
    will be executed and, on exit from the block, all warnings issued
    during execution of the block will be checked to see if the warning
    specified was issued.

    :class:`assertWarns` requires that the matched warning be issued *exactly
    once* within the context manager, or the test function will fail (whether
    the warning was issued not at all, or multiple times).
    :class:`assertDoesntWarn` requires that matched warnings are not issued
    at all.

    All other warnings are issued as normal, although the warning will not
    be shown (e.g. printed) until the exit of the context manager. If
    code within the context manager issues an exception, the test for
    warnings will be skipped (test failure will not be issued), and all
    warnings shown before the exception propagates up the stack.

    The parameters determining which warning to match are for the code
    referenced in the warning, not necessarily the code being warned.
    E.g. if code calls a deprecated function, and the deprecated function
    issues a ``DeprecationWarning``, what is matched may be the code in the
    deprecated function, not the caller; see the ``stacklevel`` parameter
    to :func:`~warnings.warn` for how this may be changed.

    Parameters
    ----------
    test : unittest.TestCase
        The test case from which this is being called, almost always
        ``self`` (so the :meth:`~unittest.TestCase.fail` method is available).

    action : {'always', ``None``, 'default', 'error', 'ignore', 'module',
              'once'}
        Unless ``None``, a warning filter matching the specified warning will
        be added to the filter before executing the block. 'always'
        (default) is generally recommended to make sure the tested
        warning will be raised. If 'always' is specified, on Python 2 the log
        of previously-issued warnings will be edited to work around a
        `Python bug <https://stackoverflow.com/questions/56821539/>`_. In this
        case using ``module`` is strongly recommended to minimize the impact
        of this editing. This filter will be removed on completion of the
        block.

    message : str, optional
        Regular expression to match the start of warning message. Default
        is to match all.

    category : type, optional
        Matches if the warning is an instance of this type or a subclass.
        Default is the base `Warning` type (matches all warnings).

    module : str, optional
        Regular expression to match the start of the name of the module
        from which the warning is issued. This is primarily used in setting
        up the warnings filter; matching the module to determine if the
        desired warning was issued is based on filename and may not work
        for modules built into the interpreter. Default is to match all.

    lineno : int, optional
        Line number from which the warning was issued. This is rarely
        useful since it will change from version to version. Default (0) will
        match all lines.

    See Also
    --------
    warnings.filterwarnings :
        The ``action``, ``message``, ``category``, ``module``, and
        ``lineno`` parameters are based on the filter specifications.

    Examples
    --------
    This class is primarily useful as part of a test suite, and cannot
    be easily demonstrated through interactive examples. See the tests
    of it in ``test_testing.py`` and its usage throughout the test suite.
    """
    requireWarning = True
    """If True, requires that the matching warning be issued (i.e. fail
       if the warning isn't issued.) If False, fail if warning is issued."""
    def __init__(self, test, action='always', message='', category=Warning,
                 module='', lineno=0):
        self._filterspec = (action, message, category, module, lineno)
        self._testcase = test
        super(assertWarns, self).__init__(record=True)

    def __enter__(self):
        """Enter the context manager, called at start of block."""
        self._log = super(assertWarns, self).__enter__()
        """Log of warnings issued within context block."""
        if self._filterspec[0] is not None:
            warnings.filterwarnings(*self._filterspec)
        if self._filterspec[0] == 'always' and sys.version_info[0:2] == (2, 7):
            # Bug in 2.7: 'always' doesn't work if warning was previously
            # issued, so remove any record of it being issued, which
            # is stored by module.
            msg_pat = re.compile(self._filterspec[1], re.I)
            cat = self._filterspec[2]
            mod_pat = re.compile(self._filterspec[3])
            for m in list(sys.modules):
                if mod_pat.match(m)\
                   and hasattr(sys.modules[m], '__warningregistry__'):
                    reg = sys.modules[m].__warningregistry__
                    for k in list(reg.keys()):
                        if msg_pat.match(k[0]) and issubclass(k[1], cat):
                            del reg[k]
                            break

    def __exit__(self, *exc_info):
        """Exit context manager, called at exit of block"""
        retval = super(assertWarns, self).__exit__(*exc_info)
        if exc_info[0] is not None: # Exception in handling, show all warnings
            for w in self._log:
                warnings.showwarning(
                    w.message, w.category, w.filename, w.lineno)
            return retval
        n_match = 0 # Number of matching warnings
        msg_pat = re.compile(self._filterspec[1], re.I)
        cat = self._filterspec[2]
        mod_pat = self._filterspec[3]
        matchall = not mod_pat #Empty pattern, match all modules
        mod_pat = re.compile(mod_pat)
        lineno = int(self._filterspec[4])
        # show_warnings isn't given the module, just the filename,
        # so find filenames of desired modules.
        mod_files = []
        for m in list(sys.modules):
            if mod_pat.match(m) and hasattr(sys.modules[m], '__file__'):
                fnl = sys.modules[m].__file__
                if fnl is None:
                    continue
                if fnl.lower().endswith(('.pyc', '.pyo')):
                    fnl = fnl[:-1]
                mod_files.append(fnl)
        for w in self._log:
            if issubclass(w.category, cat) \
               and (matchall or w.filename in mod_files) \
               and msg_pat.match(str(w.message)) and lineno in (0, w.lineno):
                n_match += 1
            else:
                warnings.showwarning(
                    w.message, w.category, w.filename, w.lineno)
        if self.requireWarning:
            if n_match == 0:
                self._testcase.fail('Warning not issued.')
            elif n_match > 1:
                self._testcase.fail('Warning issued {} times.'.format(n_match))
        elif n_match:
            self._testcase.fail('Warning was issued.')


class assertDoesntWarn(assertWarns):
    __doc__ = 'Tests that a warning is not raised.' \
              + assertWarns.__doc__[assertWarns.__doc__.index('\n'):]
    requireWarning = False
    pass


add_build_to_path()
