#!/usr/bin/env python

"""Support for ISTP-compliant CDFs

The `ISTP metadata standard  <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_
specifies the interpretation of the attributes in a CDF to describe
relationships between the variables and their physical interpretation.

This module supports that subset of CDFs.

Authors: Jon Niehof

Additional Contributors: Lorna Ellis, Asher Merrill

Institution: University of New Hampshire

Contact: Jonathan.Niehof@unh.edu

.. rubric:: Classes

.. autosummary::
    :toctree: autosummary  
    :template: clean_class.rst

    FileChecks
    VariableChecks

.. rubric:: Functions

.. autosummary::
    :toctree: autosummary  

    fillval
    format
"""

import datetime
import math
import os.path
import re

import numpy
import spacepy.datamodel
import spacepy.pycdf
import spacepy.pycdf.const


class VariableChecks(object):
    """ISTP compliance checks for a single variable.

    Checks a variable's compliance with ISTP standards. This mostly
    performs checks that are not currently performed by the `ISTP
    skeleton editor <https://spdf.gsfc.nasa.gov/skteditor/>`_.  All
    tests return a list, one error string for every noncompliance
    found (empty list if compliant). :meth:`all` will perform all
    tests and concatenate all errors.

    .. autosummary::

        all
        depends
        depsize
        fieldnam
        recordcount
        validplottype
        validrange
        validscale
        
    .. automethod:: all
    .. automethod:: depends
    .. automethod:: depsize
    .. automethod:: fieldnam
    .. automethod:: recordcount
    .. automethod:: validplottype
    .. automethod:: validrange
    .. automethod:: validscale

    """
    #When adding new tests, add to list above, and the list in all()
    #Validation failures should be formatted as a sentence (initial cap,
    #closing period) and NOT include the variable name.

    @classmethod
    def all(cls, v):
        """Perform all variable tests

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> f = spacepy.pycdf.CDF('foo.cdf', create=True)
        >>> v = f.new('Var', data=[1, 2, 3])
        >>> spacepy.pycdf.istp.VariableChecks.all(v)
        ['No FIELDNAM attribute.']
        """
        #Update this list when adding new test functions
        callme = (cls.depends, cls.depsize, cls.fieldnam, cls.recordcount,
                  cls.validrange, cls.validscale, cls.validplottype)
        errors = []
        for f in callme:
            errors.extend(f(v))
        return errors

    @classmethod
    def depends(cls, v):
        """Checks that DEPEND and LABL_PTR variables actually exist

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        return ['{} variable {} missing'.format(a, v.attrs[a])
                for a in v.attrs
                if a.startswith(('DEPEND_', 'LABL_PTR_')) and
                not v.attrs[a] in v.cdf_file]

    @classmethod
    def depsize(cls, v):
        """Checks that DEPEND has same shape as that dim

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        rv = int(v.rv()) #RV is a leading dimension
        errs = []
        # Check that don't have invalid DEPEND_1
        if v.shape == (0,):
            if 'DEPEND_1' in v.attrs or 'DEPEND_2' in v.attrs:
                errs.append('Do not expect DEPEND_1 or DEPEND_2 in 1 dimensional variable.')
        for i in range(rv, len(v.shape)): #This is index on shape (of var)
            depidx = i + 1 - rv #This is x in  DEPEND_x
            target = v.shape[i]
            if not 'DEPEND_{}'.format(depidx) in v.attrs:
                continue
            d = v.attrs['DEPEND_{}'.format(depidx)]
            if d in v.cdf_file:
                dv = v.cdf_file[d]
            else:
                continue #this is a different error
            #We hope the only weirdness is whether the dependency
            #is constant, or dependent on record. If it's dependent
            #on another dependency, this gets really weird really fast
            # If the dependency is dependent, remove the lower level
            # dependency size from consideration
            # eg. if counts [80,48], depends on energy [80,48],
            # depends on look [80], remove 80 from the view of energy
            # so that we accurately check 48==48.
            # NB: This assumes max of two layers of dependency
            if 'DEPEND_2' in dv.attrs:
                errs.append('Do not expect three layers of dependency.')
                continue
            elif 'DEPEND_1' in dv.attrs:
                dd = dv.attrs['DEPEND_1']
                if dd in v.cdf_file:
                    ddv = v.cdf_file[dd]
                else:
                    continue #this is a different error
                actual = list(dv.shape)
                for ii in actual:
                    if ii in ddv.shape:
                        actual.remove(ii)
                if 'DEPEND_0' in dv.attrs:
                    # record varying
                    dd = dv.attrs['DEPEND_0']
                    if dd[:5] != 'Epoch':
                        errs.append('Expect DEPEND_0 to be Epoch.')
                        continue
                    if dd in v.cdf_file:
                        ddv = v.cdf_file[dd]
                    else:
                        continue #this is a different error
                    for ii in actual:
                        if ii in ddv.shape:
                            actual.remove(ii)
                    
                if len(actual) != 1:
                    errs.append('More complicated double dependency than taken into account.')
                    continue
                else:
                    actual = actual[0]
            else:
                actual = dv.shape[int(dv.rv())]
            if target != actual:
                errs.append('Dim {} sized {} but DEPEND_{} {} sized {}.'.format(
                    i, target, depidx, d, actual))

        return errs

    @classmethod
    def recordcount(cls, v):
        """Check that the DEPEND_0 has same record count as variable

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        if not v.rv() or not 'DEPEND_0' in v.attrs:
            return []
        dep0 = v.attrs['DEPEND_0']
        if not dep0 in v.cdf_file: #This is a DIFFERENT error
            return []
        if len(v) != len(v.cdf_file[dep0]):
            return ['{} records; DEPEND_0 {} has {}.'.format(
                len(v), dep0, len(v.cdf_file[dep0]))]
        return []

    @classmethod
    def validrange(cls, v):
        """Check that all values are within VALIDMIN/VALIDMAX, or FILLVAL

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        #This makes datetime comparisons a lot easier
        raw_v = v.cdf_file.raw_var(v.name())
        vshape = v.shape
        minval, maxval = spacepy.pycdf.lib.get_minmax(v.type())
        data = raw_v[...]
        is_fill = numpy.isclose(data, raw_v.attrs['FILLVAL']) \
                  if 'FILLVAL' in raw_v.attrs \
                     else numpy.zeros(shape=data.shape, dtype=numpy.bool)
        for which in ('VALIDMIN', 'VALIDMAX'):
            if not which in v.attrs:
                continue
            #Some of these comparisons better raw, but display non-raw
            rawattrval = raw_v.attrs[which]
            attrval = v.attrs[which]
            multidim = bool(numpy.shape(attrval)) #multi-dimensional
            if multidim: #Compare shapes, require only 1D var
                #Match attribute dim to first non-record var dim
                firstdim = int(v.rv())
                if vshape[firstdim] != numpy.shape(attrval)[0]:
                    errs.append(('{} element count {} does not match first data'
                                 ' dimension size {}.').format(
                                     which, numpy.shape(attrval)[0],
                                     data.shape[firstdim]))
                    continue
                if len(vshape) != firstdim + 1: #only one non-record dim
                    errs.append('Multi-element {} only valid with 1D variable.'
                                .format(which))
                    continue
                if firstdim: #Add pseudo-record dim
                    attrval = numpy.reshape(attrval, (1, -1))
            if numpy.any((rawattrval < minval)) \
               or numpy.any((rawattrval > maxval)):
                errs.append('{} ({}) outside data range ({},{}) for {}.'.format(
                    which, rawattrval, minval, maxval, v.name()))
            if not len(data): #nothing to compare
                continue
            bigger, smaller = (data, rawattrval) if which == 'VALIDMIN' \
                              else (rawattrval, data)
            idx = bigger < smaller
            idx = numpy.logical_and(idx, numpy.logical_not(is_fill))
            if idx.any():
                badidx = numpy.nonzero(idx)
                badvals = v[...][badidx] #non-raw data for display
                if len(badidx) > 1: #Multi-dimensional data
                    badidx = numpy.transpose(badidx) #Group by value not axis
                else:
                    badidx = badidx[0] #Just recover the index value
                direction = 'under' if which == 'VALIDMIN' else 'over'
                errs.append('Value {} at index {} {} {} {}.'.format(
                    ', '.join(str(d) for d in badvals),
                    ', '.join(str(d) for d in badidx),
                    direction, which,
                    attrval[0, :] if multidim else attrval))
        if ('VALIDMIN' in v.attrs) and ('VALIDMAX' in v.attrs):
            if numpy.any(v.attrs['VALIDMIN'] > v.attrs['VALIDMAX']):
                errs.append('VALIDMIN > VALIDMAX.')
        return errs

    @classmethod
    def validscale(cls, v):
        """Check that SCALEMIN<=SCALEMAX, and neither goes out 
        of range for CDF datatype.

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        #This makes datetime comparisons a lot easier
        raw_v = v.cdf_file.raw_var(v.name())
        vshape = v.shape
        minval, maxval = spacepy.pycdf.lib.get_minmax(raw_v.type())
        for which in ('SCALEMIN', 'SCALEMAX'):
            if not which in v.attrs:
                continue
            rawattrval = raw_v.attrs[which]
            #This is all c/p from validrange, should pull out
            #but they're slightly different contexts.
            multidim = bool(numpy.shape(rawattrval)) #multi-dimensional
            if multidim: #Compare shapes, require only 1D var
                #Match attribute dim to first non-record var dim
                firstdim = int(v.rv())
                if vshape[firstdim] != numpy.shape(rawattrval)[0]:
                    errs.append(('{} element count {} does not match first data'
                                 ' dimension size {}.').format(
                                     which, numpy.shape(rawattrval)[0],
                                     v.shape[firstdim]))
                    continue
                if len(vshape) != firstdim + 1: #only one non-record dim
                    errs.append('Multi-element {} only valid with 1D variable.'
                                .format(which))
                    continue
                if firstdim: #Add pseudo-record dim
                    rawattrval = numpy.reshape(rawattrval, (1, -1))
            if numpy.any(rawattrval < minval) or numpy.any(rawattrval > maxval):
                errs.append('{} ({}) outside data range ({},{}).'.format(
                    which, rawattrval[0, :] if multidim else rawattrval,
                    minval, maxval))
        if ('SCALEMIN' in raw_v.attrs) and ('SCALEMAX' in raw_v.attrs):
            if numpy.any(raw_v.attrs['SCALEMIN'] > raw_v.attrs['SCALEMAX']):
                errs.append('SCALEMIN > SCALEMAX.')
        return errs

    @classmethod
    def validplottype(cls, v):
        """Check that plottype matches dimensions.

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        time_st = 'time_series'
        spec_st = 'spectrogram'
        errs = []
        if 'DISPLAY_TYPE' in v.attrs:
            if (len(v.shape) == 1) and (v.attrs['DISPLAY_TYPE'] != time_st):
                errs.append('1 dim variable with {} display type.'.format(
                    v.attrs['DISPLAY_TYPE']))
            elif (len(v.shape) > 1) and (v.attrs['DISPLAY_TYPE'] != spec_st):
                errs.append('Multi dim variable with {} display type.'.format(
                    v.attrs['DISPLAY_TYPE']))
        return errs

    @classmethod
    def fieldnam(cls, v):
        """Check that FIELDNAM attribute matches variable name.

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        vname = v.name()
        if 'FIELDNAM' not in v.attrs:
            errs.append('No FIELDNAM attribute.')
        elif v.attrs['FIELDNAM'] != vname:
            errs.append('FIELDNAM attribute {} does not match var name.'
                        .format(v.attrs['FIELDNAM']))
        return errs


class FileChecks(object):
    """ISTP compliance checks for a CDF file.

    Checks a file's compliance with ISTP standards. This mostly
    performs checks that are not currently performed by the `ISTP
    skeleton editor <https://spdf.gsfc.nasa.gov/skteditor/>`_.  All
    tests return a list, one error string for every noncompliance
    found (empty list if compliant). :meth:`all` will perform all
    tests and concatenate all errors.

    .. autosummary::

        all
        filename
        time_monoton
        times
        
    .. automethod:: all
    .. automethod:: filename
    .. automethod:: time_monoton
    .. automethod:: times

    """
    #When adding new tests, add to list above, and the list in all()
    #Validation failures should be formatted as a sentence (initial cap,
    #closing period).

    @classmethod
    def all(cls, f):
        """Perform all variable and file-level tests

        In addition to calling every test in this class, will also call
        :meth:`VariableChecks.all` for every variable in the file.

        Parameters
        ----------
        f : :class:`~spacepy.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> f = spacepy.pycdf.CDF('foo.cdf', create=True)
        >>> v = f.new('Var', data=[1, 2, 3])
        >>> spacepy.pycdf.istp.FileChecks.all(f)
        ['No Logical_source in global attrs.',
        'No Logical_file_id in global attrs.',
        'Cannot parse date from filename foo.cdf.',
        'Var: No FIELDNAM attribute.']
        """
        #Update this list when adding new test functions
        callme = (cls.filename, cls.time_monoton, cls.times,)
        errors = []
        for func in callme:
            errors.extend(func(f))
        for v in f:
            errors.extend(('{}: {}'.format(v, e)
                           for e in VariableChecks.all(f[v])))
        return errors
                
    @classmethod
    def filename(cls, f):
        """Compare filename to global attributes

        Check that the logical_file_id matches the actual filename,
        and logical_source also matches.

        Parameters
        ----------
        f : :class:`~spacepy.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        for a in ('Logical_source', 'Logical_file_id'):
            if not a in f.attrs:
                errs.append('No {} in global attrs.'.format(a))
        if errs:
            return errs
        fname = os.path.basename(f.pathname)
        if not bytes is str:
            fname = fname.decode('ascii')
        if not fname.startswith(f.attrs['Logical_source'][0]):
            errs.append("Logical_source {} doesn't match filename {}.".format(
                f.attrs['Logical_source'][0], fname))
        if fname[:-4] != f.attrs['Logical_file_id'][0]:
            errs.append("Logical_file_id {} doesn't match filename {}.".format(
                f.attrs['Logical_file_id'][0], fname))
        return errs

    @classmethod
    def time_monoton(cls, f):
        """Checks that times are monotonic

        Check that all Epoch variables are monotonically increasing.

        Parameters
        ----------
        f : :class:`~spacepy.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        for v in f:
            if not f[v].type() in (spacepy.pycdf.const.CDF_EPOCH.value,
                                   spacepy.pycdf.const.CDF_EPOCH16.value,
                                   spacepy.pycdf.const.CDF_TIME_TT2000.value):
                continue
            data = f[v][...]
            idx = numpy.where(numpy.diff(data) < datetime.timedelta(0))[0]
            if not any(idx):
                continue
            errs.append('{}: Nonmonotonic time at record {}.'.format(
                v, ', '.join((str(i) for i in (idx + 1)))))
        return errs

    @classmethod
    def times(cls, f):
        """Compare filename to times

        Check that all Epoch variables only contain times matching filename

        Parameters
        ----------
        f : :class:`~spacepy.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.

        Notes
        -----
        This function assumes daily files and should be extended based on the
        File_naming_convention global attribute (which itself is another good
        check to have.)
        """
        errs = []
        fname = os.path.basename(f.pathname)
        if not bytes is str:
            fname = fname.decode('ascii')
        m = re.search('\d{8}', fname)
        if not m:
            return ['Cannot parse date from filename {}'.format(fname)]
        datestr = m.group(0)
        for v in f:
            if f[v].type() in (spacepy.pycdf.const.CDF_EPOCH.value,
                               spacepy.pycdf.const.CDF_EPOCH16.value,
                               spacepy.pycdf.const.CDF_TIME_TT2000.value):
                datestrs = list(set((d.strftime('%Y%m%d') for d in f[v][...])))
                if len(datestrs) == 0:
                    continue
                elif len(datestrs) > 1:
                    errs.append('{}: multiple days {}.'.format(
                        v, ', '.join(sorted(datestrs))))
                elif datestrs[0] != datestr:
                    errs.append('{}: date {} doesn\'t match file {}.'.format(
                        v, datestrs[0], fname))
        return errs


def fillval(v): 
    """Set ISTP-compliant FILLVAL on a variable

    Sets a CDF variable's `FILLVAL
    <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL>`_
    attribute to the value required by ISTP (based on variable type).

    Parameters
    ----------
    v : :class:`~spacepy.pycdf.Var`
        CDF variable to update

    Examples
    --------
    >>> import spacepy.pycdf
    >>> import spacepy.pycdf.istp
    >>> f = spacepy.pycdf.CDF('foo.cdf', create=True)
    >>> v = f.new('Var', data=[1, 2, 3])
    >>> spacepy.pycdf.istp.fillval(v)
    >>> v.attrs['FILLVAL']
    -128
    """
    #Fill value, indexed by the CDF type (numeric)
    fillvals = {}
    #Integers
    for i in (1, 2, 4, 8):
        fillvals[getattr(spacepy.pycdf.const, 'CDF_INT{}'.format(i)).value] = \
            - 2 ** (8*i - 1)
        if i == 8:
            continue
        fillvals[getattr(spacepy.pycdf.const, 'CDF_UINT{}'.format(i)).value] = \
            2 ** (8*i) - 1
    fillvals[spacepy.pycdf.const.CDF_EPOCH16.value] = (-1e31, -1e31)
    fillvals[spacepy.pycdf.const.CDF_REAL8.value] = -1e31
    fillvals[spacepy.pycdf.const.CDF_REAL4.value] = -1e31
    fillvals[spacepy.pycdf.const.CDF_CHAR.value] = ' '
    fillvals[spacepy.pycdf.const.CDF_UCHAR.value] = ' '
    #Equivalent pairs
    for cdf_t, equiv in (
            (spacepy.pycdf.const.CDF_TIME_TT2000, spacepy.pycdf.const.CDF_INT8),
            (spacepy.pycdf.const.CDF_EPOCH, spacepy.pycdf.const.CDF_REAL8),
            (spacepy.pycdf.const.CDF_BYTE, spacepy.pycdf.const.CDF_INT1),
            (spacepy.pycdf.const.CDF_FLOAT, spacepy.pycdf.const.CDF_REAL4),
            (spacepy.pycdf.const.CDF_DOUBLE, spacepy.pycdf.const.CDF_REAL8),
    ):
        fillvals[cdf_t.value] = fillvals[equiv.value]
    if 'FILLVAL' in v.attrs:
        del v.attrs['FILLVAL']
    v.attrs.new('FILLVAL', data=fillvals[v.type()], type=v.type())


def format(v, use_scaleminmax=False, dryrun=False):
    """Set ISTP-compliant FORMAT on a variable

    Sets a CDF variable's `FORMAT
    <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FORMAT>`_
    attribute, which provides a Fortran-like format string that should
    be useable for printing any valid value in the variable. Sets
    according to the VALIDMIN/VALIDMAX attributes (or, optionally,
    SCALEMIN/SCALEMAX) if present, otherwise uses the full range of
    the type.

    Parameters
    ----------
    v : :class:`~spacepy.pycdf.Var`
        Variable to update
    use_scaleminmax : bool, optional
        Use SCALEMIN/MAX instead of VALIDMIN/MAX (default False).
        Note: istpchecks may complain about result.
    dryrun : bool, optional
        Print the decided format to stdout instead of modifying
        the CDF (for use in command-line debugging) (default False).

    Examples
    --------
    >>> import spacepy.pycdf
    >>> import spacepy.pycdf.istp
    >>> f = spacepy.pycdf.CDF('foo.cdf', create=True)
    >>> v = f.new('Var', data=[1, 2, 3])
    >>> spacepy.pycdf.istp.format(v)
    >>> v.attrs['FORMAT']
    'I4'

    """
    if use_scaleminmax:
        minn = 'SCALEMIN'
        maxx = 'SCALEMAX'
    else:
        minn = 'VALIDMIN'
        maxx = 'VALIDMAX'
    cdftype = v.type()
    if cdftype in (spacepy.pycdf.const.CDF_INT1.value,
                   spacepy.pycdf.const.CDF_INT2.value,
                   spacepy.pycdf.const.CDF_INT4.value,
                   spacepy.pycdf.const.CDF_INT8.value,
                   spacepy.pycdf.const.CDF_UINT1.value,
                   spacepy.pycdf.const.CDF_UINT2.value,
                   spacepy.pycdf.const.CDF_UINT4.value,
                   spacepy.pycdf.const.CDF_BYTE.value):
        if minn in v.attrs: #Just use validmin or scalemin
            minval = v.attrs[minn]
        elif cdftype in (spacepy.pycdf.const.CDF_UINT1.value,
                         spacepy.pycdf.const.CDF_UINT2.value,
                         spacepy.pycdf.const.CDF_UINT4.value): #unsigned, easy
            minval = 0
        elif cdftype == spacepy.pycdf.const.CDF_BYTE.value:
            minval = - 2 ** 7
        else: #Signed, harder
            size = next((i for i in (1, 2, 4, 8) if getattr(
                spacepy.pycdf.const, 'CDF_INT{}'.format(i)).value == cdftype))
            minval = - 2 ** (8*size  - 1)
        if maxx in v.attrs: #Just use max
            maxval = v.attrs[maxx]
        elif cdftype == spacepy.pycdf.const.CDF_BYTE.value:
            maxval = 2 ** 7 - 1
        else:
            size = next((8 * i for i in (1, 2, 4) if getattr(
                spacepy.pycdf.const, 'CDF_UINT{}'.format(i)).value == cdftype),
                        None)
            if size is None:
                size = next((8 * i for i in (1, 2, 4, 8) if getattr(
                    spacepy.pycdf.const, 'CDF_INT{}'.format(i)).value ==
                             cdftype)) - 1
            maxval = 2 ** size - 1
        #Two tricks:
        #-Truncate and add 1 rather than ceil so get
        #powers of 10 (log10(10) = 1 but needs two digits)
        #-Make sure not taking log of zero
        if minval < 0: #Need an extra space for the negative sign
            fmt = 'I{}'.format(int(math.log10(max(
                abs(maxval), abs(minval), 1))) + 2)
        else:
            fmt = 'I{}'.format(int(
                math.log10(maxval) if maxval != 0 else 1) + 1)
    elif cdftype == spacepy.pycdf.const.CDF_TIME_TT2000.value:
        fmt = 'A{}'.format(len('9999-12-31T23:59:59.999999999'))
    elif cdftype == spacepy.pycdf.const.CDF_EPOCH16.value:
        fmt = 'A{}'.format(len('31-Dec-9999 23:59:59.999.999.000.000'))
    elif cdftype == spacepy.pycdf.const.CDF_EPOCH.value:
        fmt = 'A{}'.format(len('31-Dec-9999 23:59:59.999'))
    elif cdftype in (spacepy.pycdf.const.CDF_REAL8.value,
                     spacepy.pycdf.const.CDF_REAL4.value,
                     spacepy.pycdf.const.CDF_FLOAT.value,
                     spacepy.pycdf.const.CDF_DOUBLE.value):
        # Prioritize SCALEMIN/MAX to find the number of decimals to include
        if 'SCALEMIN' in v.attrs and 'SCALEMAX' in v.attrs:
            range = v.attrs['SCALEMAX'] - v.attrs['SCALEMIN']
        # If not, use VALIDMIN/MAX
        elif 'VALIDMIN' in v.attrs and 'VALIDMAX' in v.attrs:
            range = v.attrs['VALIDMAX'] - v.attrs['VALIDMIN']
        # If not, just use nothing.
        else:
            range = None
        # Find how many spaces we need for the 'integer' part of the number
        # (Use maxx-minn for this...effectively uses VALIDMIN/MAX for most
        # cases.)
        if range and (minn in v.attrs and maxx in v.attrs):
            if len(str(int(v.attrs[maxx]))) >=\
               len(str(int(v.attrs[minn]))):
                ln = str(int(v.attrs[maxx]))
            else:
                ln = str(int(v.attrs[minn]))
        if range and ln and range < 0: # Cover all our bases:
            # raise ValueError('Range ({} - {}) cannot be negative:'
                # '\nVarname: {}\nRange: {}'.format(maxx, minn, v, range))
            ### Instead of throwing an error, just use None
            # There are old cases that for some reason have negative ranges, so
            # this is really more of a compatibility choice than a good
            # decision.
            range = None
        # All of the lengths below (the +4, +3, +2, etc...) should be EXACTLY
        # enough.  Consider adding 1, (4+1=5, 3+1=4, etc...) to possibly make
        # this easier.
        # elif range and ln and range <= 11: # If range <= 11, we want 2 decimal places:
        if range and ln and range <= 11: # If range <= 11, we want 2 decimal places:
            # Need extra for '.', and 3 decimal places (4 extra)
            fmt = 'F{}.3'.format(len([i for i in ln]) + 4)
        elif range and ln and 11 < range <= 101:
            # Need extra for '.' (1 extra)
            fmt = 'F{}.2'.format(len([i for i in ln]) + 3)
        elif range and ln and 101 < range <= 1000:
            # Need extra for '.' (1 extra)
            fmt = 'F{}.1'.format(len([i for i in ln]) + 2)
        else:
            # No range, must not be populated, copied from REAL4/8(s) above
            # OR we don't care because it's a 'big' number:
            fmt = 'G10.2E3'
    elif cdftype in (spacepy.pycdf.const.CDF_CHAR.value,
                     spacepy.pycdf.const.CDF_UCHAR.value):
        #This is a bit weird but pycdf consider nelems private. Should change.
        fmt = 'A{}'.format(v._nelems())
    else:
        raise ValueError("Couldn't find FORMAT for {} of type {}".format(
            v.name(),
            spacepy.pycdf.lib.cdftypenames.get(cdftype, 'UNKNOWN')))
    if dryrun:
        print(fmt)
    else:
        if 'FORMAT' in v.attrs:
            del v.attrs['FORMAT']
        v.attrs.new('FORMAT', data=fmt, type=spacepy.pycdf.const.CDF_CHAR)


