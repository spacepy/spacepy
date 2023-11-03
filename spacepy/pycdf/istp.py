#!/usr/bin/env python

"""Support for ISTP-compliant CDFs

The `ISTP metadata standard <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_
specifies the interpretation of the attributes in a CDF to describe
relationships between the variables and their physical interpretation.

This module supports that subset of CDFs.

Authors: Jon Niehof

Additional Contributors: Lorna Ellis, Asher Merrill

Institution: University of New Hampshire

Contact: Jonathan.Niehof@unh.edu

.. rubric:: Classes

.. autosummary::
    :toctree:
    :template: clean_class.rst

    FileChecks
    VarBundle
    VariableChecks

.. rubric:: Functions

.. autosummary::
    :toctree:

    fillval
    format
    nanfill
"""

import collections
import datetime
import functools
import inspect
import itertools
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
    found (empty list if compliant). `all` will perform all
    tests and concatenate all errors.

    .. autosummary::

        all
        deltas
        depends
        depsize
        empty_entry
        fieldnam
        fillval
        recordcount
        validdisplaytype
        validrange
        validscale
        
    .. automethod:: all
    .. automethod:: deltas
    .. automethod:: depends
    .. automethod:: depsize
    .. automethod:: empty_entry
    .. automethod:: fieldnam
    .. automethod:: fillval
    .. automethod:: recordcount
    .. automethod:: validdisplaytype
    .. automethod:: validrange
    .. automethod:: validscale

    """
    #When adding new tests, add to list above
    #Validation failures should be formatted as a sentence (initial cap,
    #closing period) and NOT include the variable name.

    @classmethod
    def all(cls, v, catch=False):
        """Perform all variable tests

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check
        catch : bool
            Catch exceptions in tests (default False). If True, any
            exceptions in subtests will result in an addition to the
            validation failures of the form "Test x did not complete."
            Calling the individual test will reveal the full traceback.

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
        callme = [func for name, func in inspect.getmembers(cls)
                  if not name.startswith('_') and not name.endswith('_')
                  and callable(func) and name != 'all']
        errors = []
        for f in callme:
            try:
                errors.extend(f(v))
            except:
                if catch:
                    errors.append('Test {} did not complete.'.format(
                        f.__name__))
                else:
                    raise
        return errors

    @classmethod
    def depends(cls, v):
        """Checks that DELTA, DEPEND, and LABL_PTR variables exist

        Check that variables specified in the variable attributes for
        `DELTA
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DELTA>`_,
        `DEPEND
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DEPEND_0>`_,
        and `LABL_PTR
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#LABL_PTR_1>`_
        exist in the CDF.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        return ['{} variable {} missing.'.format(a, v.attrs[a])
                for a in v.attrs
                if (a.startswith(('DEPEND_', 'LABL_PTR_',))
                    or a in ('DELTA_PLUS_VAR', 'DELTA_MINUS_VAR'))
                and not v.attrs[a] in v.cdf_file]

    @classmethod
    def deltas(cls, v):
        """Check DELTA variables

        Check that variables specified in the variable attributes for
        `DELTA
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DELTA>`_
        match the type, size, and units of this variable.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        errs = []
        if v.rv():
            shape = v.shape[1:]
            n_recs = len(v)
        else:
            shape = v.shape
            n_recs = None
        for delta in ('DELTA_PLUS_VAR', 'DELTA_MINUS_VAR'):
            if not delta in v.attrs:
                continue
            deltavar = v.cdf_file[v.attrs[delta]]
            if deltavar.type() != v.type():
                errs.append(
                    '{} type {} does not match variable type {}.'.format(
                        delta, spacepy.pycdf.lib.cdftypenames[deltavar.type()],
                        spacepy.pycdf.lib.cdftypenames[v.type()]))
            if deltavar.attrs.get('UNITS', None) != v.attrs.get('UNITS', None):
                errs.append('{} units do not match variable units.'.format(
                    delta))
            if deltavar.rv():
                dshape = deltavar.shape[1:]
                d_n_recs = len(deltavar)
            else:
                dshape = deltavar.shape
                d_n_recs = None
            if dshape != shape:
                errs.append(
                    '{} shape {} does not match variable shape {}.'.format(
                        delta, dshape, shape))
            if d_n_recs is not None and n_recs is not None \
               and d_n_recs != n_recs:
                errs.append((
                    '{} record count {} does not match variable record'
                    ' count {}.').format(
                        delta, d_n_recs, n_recs))
        return errs

    @classmethod
    def depsize(cls, v):
        """Checks that DEPEND has same shape as that dim

        Compares the size of variables specified in the variable
        attributes for `DEPEND
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DEPEND_0>`_
        and compares to the size of the corresponding dimension in
        this variable.

        Parameters
        ----------
        v : `~.pycdf.Var`
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
            if dv.rv() != ('DEPEND_0' in dv.attrs):
                errs.append('DEPEND_{} {} is RV but has no DEPEND_0.'
                            .format(depidx, d))
                continue
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
    def empty_entry(cls, v):
        """Check for attributes with empty string

        Checks attributes for this variable for any entries consisting
        of an empty string. These should be replaced with a single space.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        for a in v.attrs:
            if v.attrs.type(a) in (spacepy.pycdf.const.CDF_CHAR.value,
                                   spacepy.pycdf.const.CDF_UCHAR.value) \
                and v.attrs[a] == '':
                errs.append('Empty CHAR entry for attribute {}.'.format(a))
        return errs

    @classmethod
    def fillval(cls, v):
        """Check for FILLVAL presence, type, value

        Checks variable for existence of `FILLVAL
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL>`_
        attribute and makes sure it is the same type as variable and matches
        ISTP value.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        See Also
        --------
        spacepy.pycdf.istp.fillval : Automatic setting of this value.
        """
        errs = []
        if not 'FILLVAL' in v.attrs:
            return ['No FILLVAL attribute.']
        if v.attrs.type('FILLVAL') != v.type():
            errs.append(
                'FILLVAL type {} does not match variable type {}.'.format(
                    spacepy.pycdf.lib.cdftypenames[v.attrs.type('FILLVAL')],
                    spacepy.pycdf.lib.cdftypenames[v.type()]))
        expected = fillval(v, ret=True)
        timetype = v.type() in spacepy.pycdf.lib.timetypes
        actual = (v.cdf_file.raw_var(v.name()) if timetype else v)\
                 .attrs['FILLVAL']
        match = numpy.isclose(
            actual, expected, atol=0, rtol=1e-7)\
            if numpy.issubdtype(v.dtype, numpy.floating)\
            else numpy.all(actual == expected)
        if not match:
            if timetype:
                converted_expected = {
                    spacepy.pycdf.const.CDF_EPOCH.value:
                    spacepy.pycdf.lib.v_epoch_to_datetime,
                    spacepy.pycdf.const.CDF_EPOCH16.value:
                    spacepy.pycdf.lib.v_epoch16_to_datetime,
                    spacepy.pycdf.const.CDF_TIME_TT2000.value:
                    spacepy.pycdf.lib.v_tt2000_to_datetime
                }[v.type()](expected)
                errs.append(
                    'FILLVAL {} ({}), should be {} ({}) for variable type {}.'
                    .format(
                        actual,
                        v.attrs['FILLVAL'],
                        expected,
                        converted_expected,
                        spacepy.pycdf.lib.cdftypenames[v.type()]))
            else:
                errs.append(
                    'FILLVAL {}, should be {} for variable type {}.'.format(
                        actual, expected,
                        spacepy.pycdf.lib.cdftypenames[v.type()]))
        return errs

    @classmethod
    def recordcount(cls, v):
        """Check that the DEPEND_0 has same record count as variable

        Checks the record count of the variable specified in the
        variable attribute for `DEPEND_0
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DEPEND_0>`_
        and compares to the record count for this variable.

        Parameters
        ----------
        v : `~.pycdf.Var`
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
    def _validhelper(cls, v, rng=True):
        """Helper function for checking SCALEMIN/MAX, VALIDMIN/MAX

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        rng : bool
            Do range check (True, default) or scale check (False)

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        validscale = 'VALID' if rng else 'SCALE'
        whichmin, whichmax = ('VALIDMIN', 'VALIDMAX') if rng \
                             else ('SCALEMIN', 'SCALEMAX')
        errs = []
        vshape = v.shape
        minval, maxval = spacepy.pycdf.lib.get_minmax(v.type())
        if rng:
            data = v[...]
            is_fill = False
            if 'FILLVAL' in v.attrs:
                filldtype = spacepy.pycdf.lib.numpytypedict.get(
                    v.attrs.type('FILLVAL'), object)
                if numpy.issubdtype(v.dtype, numpy.floating) \
                   and numpy.issubdtype(filldtype, numpy.floating):
                    is_fill = numpy.isclose(data, v.attrs['FILLVAL'])
                elif numpy.can_cast(numpy.asanyarray(v.attrs['FILLVAL']),
                                    v.dtype):
                    is_fill = data == v.attrs['FILLVAL']
        for which in (whichmin, whichmax):
            if not which in v.attrs:
                continue
            if v.attrs.type(which) != v.type():
                errs.append(
                    '{} type {} does not match variable type {}.'.format(
                        which,
                        spacepy.pycdf.lib.cdftypenames[v.attrs.type(which)],
                        spacepy.pycdf.lib.cdftypenames[v.type()]))
            attrval = v.attrs[which]
            multidim = bool(numpy.shape(attrval)) #multi-dimensional
            if multidim: #Compare shapes, require only 1D var
                #Match attribute dim to first non-record var dim
                firstdim = int(v.rv())
                if vshape[firstdim] != numpy.shape(attrval)[0]:
                    errs.append(('{} element count {} does not match first data'
                                 ' dimension size {}.').format(
                                     which, numpy.shape(attrval)[0],
                                     v.shape[firstdim]))
                    continue
                if len(vshape) != firstdim + 1: #only one non-record dim
                    errs.append('Multi-element {} only valid with 1D variable.'
                                .format(which))
                    continue
                if firstdim: #Add pseudo-record dim
                    attrval = numpy.reshape(attrval, (1, -1))
            # min, max, variable data all same dtype
            if not numpy.can_cast(numpy.asanyarray(attrval),
                                  numpy.asanyarray(minval).dtype):
                errs.append(
                    '{} type {} not comparable to variable type {}.'.format(
                        which,
                        spacepy.pycdf.lib.cdftypenames[v.attrs.type(which)],
                        spacepy.pycdf.lib.cdftypenames[v.type()]
                    ))
                continue # Cannot do comparisons
            if numpy.any((minval > attrval)) or numpy.any((maxval < attrval)):
                errs.append('{} ({}) outside valid data range ({},{}).'.format(
                    which, attrval[0, :] if multidim else attrval,
                    minval, maxval))
            if not rng or not len(v): #nothing to compare
                continue
            #Always put numpy array on the left so knows to do element compare
            idx = (data < attrval) if which == whichmin \
                  else (data > attrval)
            idx = numpy.logical_and(idx, numpy.logical_not(is_fill))
            if idx.any():
                direction = 'under' if which == whichmin else 'over'
                if len(vshape) == 0: #Scalar
                    errs.append('Value {} {} {} {}.'.format(
                        data, direction, which,
                        attrval[0, :] if multidim else attrval))
                    continue
                badidx = numpy.nonzero(idx)
                badvals = data[badidx]
                if len(badidx) > 1: #Multi-dimensional data
                    badidx = numpy.transpose(badidx) #Group by value not axis
                else:
                    badidx = badidx[0] #Just recover the index value
                if len(badvals) < 10:
                    badvalstr = ', '.join(str(d) for d in badvals)
                    badidxstr = ', '.join(str(d) for d in badidx)
                    errs.append('Value {} at index {} {} {} {}.'.format(
                        badvalstr, badidxstr,
                        direction, which,
                        attrval[0, :] if multidim else attrval))
                else:
                    errs.append('{} values {} {} {}'.format(
                        len(badvals), direction, which,
                        attrval[0, :] if multidim else attrval))
        if (whichmin in v.attrs) and (whichmax in v.attrs):
            if numpy.any(v.attrs[whichmin] > v.attrs[whichmax]):
                errs.append('{} > {}.'.format(whichmin, whichmax))
        return errs

    @classmethod
    def validrange(cls, v):
        """Check that all values are within VALIDMIN/VALIDMAX, or FILLVAL

        Compare all values of this variable to `VALIDMIN
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#VALIDMIN>`_
        and ``VALIDMAX``; fails validation if any values are below
        VALIDMIN or above ``VALIDMAX`` unless equal to `FILLVAL
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL>`_.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        return cls._validhelper(v)

    @classmethod
    def validscale(cls, v):
        """Check SCALEMIN<=SCALEMAX, and both in range for CDF datatype.

        Compares `SCALEMIN
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#SCALEMIN>`_
        to ``SCALEMAX`` to make sure it isn't larger and both are
        within range of the variable CDF datatype.

        Parameters
        ----------
        v : `~.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        return cls._validhelper(v, False)

    @classmethod
    def validdisplaytype(cls, v):
        """Check that plottype matches dimensions.

        Check `DISPLAYTYPE
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#DISPLAY_TYPE>`_
        of this variable and makes sure it is reasonable for the
        variable dimensions.

        Parameters
        ----------
        v : `~.pycdf.Var`
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

        Compare `FIELDNAM
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FIELDNAM>`_
        attribute to the variable name; fail validation if they don't
        match.

        Parameters
        ----------
        v : `~.pycdf.Var`
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
    found (empty list if compliant). `all` will perform all
    tests and concatenate all errors.

    .. autosummary::

        all
        empty_entry
        filename
        time_monoton
        times
        
    .. automethod:: all
    .. automethod:: empty_entry
    .. automethod:: filename
    .. automethod:: time_monoton
    .. automethod:: times

    """
    #When adding new tests, add to list above.
    #Validation failures should be formatted as a sentence (initial cap,
    #closing period).

    @classmethod
    def all(cls, f, catch=False):
        """Perform all variable and file-level tests

        In addition to calling every test in this class, will also call
        `VariableChecks.all` for every variable in the file.

        Parameters
        ----------
        f : `~.pycdf.CDF`
            Open CDF file to check
        catch : bool
            Catch exceptions in tests (default False). If True, any
            exceptions in subtests will result in an addition to the
            validation failures of the form "Test x did not complete."
            Calling the individual test will reveal the full traceback.

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
        callme = [func for name, func in inspect.getmembers(cls)
                  if not name.startswith('_') and not name.endswith('_')
                  and callable(func) and name != 'all']
        errors = []
        for func in callme:
            try:
                errors.extend(func(f))
            except:
                if catch:
                    errors.append('Test {} did not complete.'.format(
                        func.__name__))
                else:
                    raise

        for v in f:
            errors.extend(('{}: {}'.format(v, e)
                           for e in VariableChecks.all(f[v], catch=catch)))
        return errors
                
    @classmethod
    def empty_entry(cls, f):
        """Check for attributes with empty string

        Checks global attributes for this variable for any entries consisting
        of an empty string. These should be replaced with a single space.

        Parameters
        ----------
        f : `~.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        errs = []
        for a in f.attrs:
            attr = f.attrs[a]
            for i in range(attr.max_idx() + 1):
                if attr.has_entry(i) \
                    and attr.type(i) in (spacepy.pycdf.const.CDF_CHAR.value,
                                         spacepy.pycdf.const.CDF_UCHAR.value) \
                    and attr[i] == '':
                    errs.append('Empty CHAR entry {} for attribute {}.'
                                .format(i, a))
        return errs

    @classmethod
    def filename(cls, f):
        """Compare filename to global attributes

        Check global attribute `Logical_file_id
        <https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html#Logical_file_id>`_
        and `Logical_source
        <https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html#Logical_source>`_
        for consistency with CDF filename.

        Parameters
        ----------
        f : `~.pycdf.CDF`
            Open CDF file to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        errs = []
        for a in ('Logical_source', 'Logical_file_id'):
            if not a in f.attrs or len(f.attrs[a]) == 0:
                errs.append('No {} in global attrs.'.format(a))
        if errs:
            return errs
        fname = os.path.basename(f.pathname)
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

        Check that all `Epoch
        <https://spdf.gsfc.nasa.gov/istp_guide/variables.html#support_data_eg1>`_
        variables are monotonically increasing.

        Parameters
        ----------
        f : `~.pycdf.CDF`
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

        Check that all `Epoch
        <https://spdf.gsfc.nasa.gov/istp_guide/variables.html#support_data_eg1>`_
        variables only contain times matching filename.

        Parameters
        ----------
        f : `~.pycdf.CDF`
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
        fname = fname.decode('ascii')
        m = re.search(r'\d{8}', fname)
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


def fillval(v, ret=False):
    """Set ISTP-compliant FILLVAL on a variable

    Sets or returns a CDF variable's `FILLVAL
    <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL>`_
    attribute to the value required by ISTP (based on variable type).

    Parameters
    ----------
    v : `~.pycdf.Var`
        CDF variable to update

    Other Parameters
    ----------------
    ret : boolean
        If True, return the value instead of setting it (Default False, set).

    Returns
    -------
    various
        If ``ret`` is True, returns the correct value for variable type (which
        may be of various Python types).  Otherwise sets the value and returns
        ``None``.

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
    value = fillvals[v.type()]
    if ret:
        return value
    if 'FILLVAL' in v.attrs:
        del v.attrs['FILLVAL']
    v.attrs.new('FILLVAL', data=value, type=v.type())


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
    v : `~.pycdf.Var`
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
        fmt = 'A{}'.format(v.nelems())
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


def nanfill(v):
    """Set fill values to NaN

    Finds all values which are equal to ``FILLVAL``, greater than
    ``VALIDMAX``, or less than ``VALIDMIN``, and replace with ``NaN``
    (not-a-number). This is an update-in-place operation; does not return
    a copy.

    Assumes a single value for ``VALIDMIN``, ``VALIDMAX``, ``FILLVAL``
    (although if the attribute is not present, will simply assume no
    restriction.)

    Only applicable to floating-point types. Best applied to a
    `~.pycdf.VarCopy` or `~.datamodel.dmarray`
    rather than `~.pycdf.Var`. Updating a variable in a CDF
    requires one write per changed value, and also will result in a CDF
    that is no longer ISTP compliant.

    Because of floating-point comparison, the matching to ``FILLVAL`` may
    fail.

    Parameters
    ----------
    v : `~.pycdf.Var` or `~.datamodel.dmarray`
        CDF variable, data, or copy to update

    Examples
    --------
    >>> import spacepy.pycdf
    >>> import spacepy.pycdf.istp
    >>> f = spacepy.pycdf.CDF('foo.cdf', create=True)
    >>> v = f.new('Var', data=[1, 2, 3, -1e31])
    >>> spacepy.pycdf.istp.fillval(v)
    >>> data = v.copy()
    >>> data
    VarCopy([1., 2., 3., -1.e31], dtype=float32)
    >>> spacepy.pycdf.istp.nanfill(data)
    >>> data
    VarCopy([1., 2., 3., nan], dtype=float32)
    """
    #If input is a zVar, read all the data; if not, this is a no-copy operation
    indata = v[...]
    badidx = numpy.zeros(shape=v.shape, dtype=bool)
    if 'FILLVAL' in v.attrs:
        badidx |= (indata == v.attrs['FILLVAL'][...])
    if 'VALIDMIN' in v.attrs:
        badidx |= (indata < v.attrs['VALIDMIN'][...])
    if 'VALIDMAX' in v.attrs:
        badidx |= (indata > v.attrs['VALIDMAX'][...])
    #Try a simple assignment with fancy indexing
    try:
        v[badidx] = numpy.nan
    except (IndexError, ValueError):
        pass
    else:
        return #success
    #Fancy indexing failed, do element-by-element assignment
    badidx = numpy.transpose(badidx.nonzero())
    for i in badidx:
        v[tuple(i)] = numpy.nan


class VarBundle(object):
    """Collective handling of ISTP-compliant variable and its dependencies.

    Representation of an ISTP-compliant variable bundled together
    with its dependencies to enable aggregate operations. Normally
    used to copy a subset of data from one CDF or SpaceData to another by
    chaining operations, or to load just the relevant data from a CDF
    into a `~.datamodel.SpaceData`.

    ``VarBundle`` operates on a single variable within a file or SpaceData
    and its various dependencies, uncertainties, labels, etc. That variable
    can be specified one of two ways. An open CDF file or
    SpaceData can be passed as the first parameter, and the name of a
    variable within it as the second parameter. Or, for CDF files, a
    :class:`~.pycdf.Var` can be passed as the only parameter, implicitly
    defining the input file (the CDF containing that variable).

    Unusual or indecipherable error messages may indicate an ISTP
    compliance issue; see `VariableChecks` for some checks.

    Parameters
    ----------
    source : `~.pycdf.CDF`, `~.datamodel.SpaceData`, or `~.pycdf.Var`
        SpaceData or open CDF containing the variable to process, or the CDF variable itself.
    name : `str`
        Name of the variable within ``source`` to process ("main variable").

    See Also
    --------
    .datamodel.fromCDF
    .pycdf.CDF.copy

    Notes
    -----
    If using :class:`~.datamodel.SpaceData` input, the contents are
    assumed to be `ISTP compliant
    <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_. In particular,
    the following attributes of the enclosed
    :class:`~.datamodel.dmarray` are used (*italics* denotes required):

        * *DEPEND_0*, *DEPEND_1*, etc.
        * LABL_PTR_0, LABL_PTR_1, etc.
        * DELTA_PLUS_VAR, DELTA_MINUS_VAR
        * VALIDMIN, VALIDMAX, *FILLVAL*

    Examples
    --------
    >>> import spacepy.pycdf
    >>> import spacepy.pycdf.istp
    >>> #https://rbsp-ect.newmexicoconsortium.org/data_pub/rbspa/hope/level3/pitchangle/2012/
    >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
    >>> infile['FPDU']
    <Var:
    CDF_FLOAT [3228, 11, 72]
    >
    >>> infile['FPDU'].attrs
    <zAttrList:
    CATDESC: HOPE differential proton flux [CDF_CHAR]
    DEPEND_0: Epoch_Ion [CDF_CHAR]
    DEPEND_1: PITCH_ANGLE [CDF_CHAR]
    DEPEND_2: HOPE_ENERGY_Ion [CDF_CHAR]
    ...
    >
    >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
    >>> b = spacepy.pycdf.istp.VarBundle(infile, 'FPDU')  # Equivalent
    >>> outfile = spacepy.pycdf.CDF('output.cdf', create=True)
    >>> b.slice(1, 2, single=True).output(outfile)
    <VarBundle:
    FPDU: CDF_FLOAT [3228, 72]
    Epoch_Ion: CDF_EPOCH [3228]
        Epoch_Ion_DELTA: CDF_REAL4 [3228]
    PITCH_ANGLE: CDF_FLOAT ---
        Pitch_LABL: CDF_CHAR*5 ---
    HOPE_ENERGY_Ion: CDF_FLOAT [3228, 72]
        ENERGY_Ion_DELTA: CDF_FLOAT [3228, 72]
        Energy_LABL: CDF_CHAR*3 [72] NRV
    >
    >>> outfile['FPDU']
    <Var:
    CDF_FLOAT [3228, 72]
    >
    >>> outfile['FPDU'].attrs
    <zAttrList:
    CATDESC: HOPE differential proton flux [CDF_CHAR]
    DEPEND_0: Epoch_Ion [CDF_CHAR]
    DEPEND_1: HOPE_ENERGY_Ion [CDF_CHAR]
    ...
    >
    >>> outfile.close()
    >>> infile.close()

    .. autosummary::

        mean
        operations
        output
        slice
        sum
        toSpaceData
        variables

    .. automethod:: mean
    .. automethod:: output
    .. automethod:: operations
    .. automethod:: slice
    .. automethod:: sum
    .. automethod:: toSpaceData
    .. automethod:: variables

    """

    def __init__(self, source, name=None):
        """Initialize variable bundle

        Parameters
        ----------
        source : `~spacepy.pycdf.CDF` or `~spacepy.pycdf.Var`
            CDF containing the variable to process, or the variable itself.
        name : `str`
            Name of the variable within ``source`` to process ("main variable").
        """
        if name is None and not hasattr(source, 'cdf_file'):
            raise TypeError('Single-argument form must be a variable'
                ' in an open CDF, not {}.'.format(type(source).__name__))
        self.mainvar = source if name is None else source[name]
        """The variable to operate on."""
        self.cdf = self.mainvar.cdf_file if name is None else source
        """Input CDF file containing the main variable."""
        self._name = self.mainvar.name() if name is None else name
        """Name of the main variable"""
        self._varinfo = {}
        """Keyed by variable name. Values are also dicts, keys are
        ``dims``, list of the main variable dimensions corresponding
        to each dimension of the variable, ``slice``, the slice
        to apply when reading this variable from the input, ``postidx``,
        a numpy fancy index to apply after reading, ``thisdim``,
        the main dimension for which this var is a dep
        (and thus it should be removed if the dim is removed),
        ``vartype``, whether this variable is the main var (M),
        a dependency (D), or DELTA of the main (U, for uncertainty),
        ``sortorder``, the order in which it should be displayed (0 for
        the main variable, 1 for dependencies, 2 for all DELTAs, and 3 for
        labels).
        """
        self._degenerate = []
        """Index by dim, is it degenerate, i.e. removed in a slice."""
        self._summed = []
        """Index by dim, is this dim summed."""
        self._mean = []
        """Index by dim, is this dim averaged."""
        self._getvarinfo()

    def _process_delta(self, mainname, deltaname):
        """Handle DELTA_PLUS/DELTA_MINUS attributes

        A DELTA variable should be the same shape and the same
        dependencies as its referrer (except potentially NRV).

        Parameters
        ----------
        mainname : str
            Name of variable that references the DELTA, i.e. it has a
            DELTA_PLUS_VAR/DELTA_MINUS_VAR attribute that references
            ``deltaname``.

        deltaname : str
            Name of the DELTA variable itself.

        Returns
        -------
        dict
            dims/slice information suitable for inclusion in ``_varinfo``.
        """
        thisvar = self.cdf[deltaname]
        mainvar = self.cdf[mainname]
        for a in thisvar.attrs: #Check that all dependencies match
            if not a.startswith(('DEPEND_', 'LABL_PTR_')):
                continue
            if a in mainvar.attrs:
                if thisvar.attrs[a] != mainvar.attrs[a]:
                    raise ValueError('{}: attribute {} mismatch with main var'
                                     .format(deltaname, a))
            elif thisvar.attrs[a] != mainname:
                raise ValueError('{}: attribute {} not in main var'
                                 .format(deltaname, a))
        rv = thisvar.rv() if hasattr(thisvar, 'rv')\
             else 'DEPEND_0' in thisvar.attrs
        if rv and not self._varinfo[mainname]['rv']:
            raise ValueError(
                '{}: Cannot handle RV DELTA with NRV variable.'
                .format(deltaname))
        thisshape = thisvar.shape
        mainshape = mainvar.shape
        if not rv and self._varinfo[mainname]['rv']: #Ignore record dim
            mainshape = mainshape[1:]
        if thisshape != mainshape:
            raise ValueError('{}: DELTA/main var shape mismatch.'
                             .format(deltaname))
        #If this is NRV and main is RV, that's okay, the R dim will
        #get removed when actually slicing.
        result = { k: self._varinfo[mainname][k][:]
                  for k in ('dims', 'slice', 'postidx') }
        result.update({
            'dv': thisvar.dv() if hasattr(thisvar, 'dv')\
                  else [True] * (len(result['dims']) - 1),
            'rv': rv,
            'sortorder': 2,
            })
        return result

    def _getvarinfo(self):
        """Find dependency and dimension information

        For main variable and its dependencies, find how dimensions
        relate to the main variable, and find all DELTA variables.
        """
        rv = self.mainvar.rv() if hasattr(self.mainvar, 'rv')\
             else 'DEPEND_0' in self.mainvar.attrs
        #Every dim maps back to itself for the main variable
        dims = list(range(len(self.mainvar.shape) + int(not rv)))
        self._degenerate = [False] * len(self.mainvar.shape)
        self._summed = [False] * len(self.mainvar.shape)
        self._mean = [False] * len(self.mainvar.shape)
        if not rv: #Fake the 0-dim
            self._degenerate.insert(0, False)
            self._summed.insert(0, False)
            self._mean.insert(0, False)
        #And every dimension is a full slice, to start
        self._varinfo[self._name] = {
            'dims': dims,
            # Dim variance is CDF concept--if not specified, assume True
            'dv': self.mainvar.dv() if hasattr(self.mainvar, 'dv')\
                  else [True] * (len(dims) - 1),
            'slice': [slice(None)] * len(dims),
            'postidx': [slice(None)] * len(dims),
            'rv': rv,
            'sortorder': 0,
            'vartype': 'M',
        }
        mainattrs = self.mainvar.attrs
        #Get the attributes that matter in the MAIN var
        attrs = {a: mainattrs[a] for a in mainattrs
                if a.startswith(('DEPEND_', 'LABL_PTR'))
                or a in ('DELTA_PLUS_VAR', 'DELTA_MINUS_VAR')}
        for a in attrs: #Process DEPEND/LABL_PTR variables
            if not a.startswith(('DEPEND_', 'LABL_PTR_')):
                continue
            thisname = attrs[a]
            if thisname in self._varinfo: #Already handled
                if self._varinfo[thisname]['sortorder'] == 3 \
                   and a.startswith('DEPEND_'):
                    #Processed before as a LABL, but also is a DEPEND.
                    #Technically ISTP violation, but have the DEPEND take
                    #priority
                    self._varinfo[thisname]['sortorder'] = 1
                continue
            thisvar = self.cdf[thisname]
            #Dimension of main var that corresponds to this var
            dim = int(a.split('_')[-1])
            dims = [0,] #Record dim always matches
            rv = thisvar.rv() if hasattr(thisvar, 'rv')\
                 else 'DEPEND_0' in thisvar.attrs or a == 'DEPEND_0'
            #For every CDF (non-record) dim, match to the main variable
            for i in range(1, len(thisvar.shape) + int(not rv)):
                #DEPEND; LABL_PTR for this dimension
                dim_dep = 'DEPEND_{}'.format(i)
                labl_dep = 'LABL_PTR_{}'.format(i)
                if not dim_dep in thisvar.attrs:
                    #No depend on this dim, so it's the dim that's represented
                    #in this variable
                    dims.append(dim)
                else: #Match to parent var
                    dim_dep = thisvar.attrs[dim_dep]
                    parentdim = next((
                        int(d.split('_')[-1]) for d in attrs
                        if (d.startswith('DEPEND_') and attrs[d] == dim_dep)
                        or (d.startswith('LABL_PTR_') and attrs[d] == labl_dep)
                    ), None)
                    if parentdim is None:
                        raise ValueError('Cannot match dim {} of {}'.format(
                            i, thisname))
                    dims.append(parentdim)
            if dims.count(dim) != 1:
                raise ValueError('Cannot find unique dimension for {}'
                                 .format(thisname))
            self._varinfo[thisname] = {
                'dims': dims,
                'dv': thisvar.dv() if hasattr(thisvar, 'dv')\
                      else [True] * (len(dims) - 1),
                'slice': [slice(None)] * len(dims),
                'postidx': [slice(None)] * len(dims),
                'rv': rv,
                'sortorder': 1 if a.startswith('DEPEND_') else 3,
                'thisdim': dim,
                'vartype': 'D',
            }
            #Process DELTAs of the DEPENDs
            for d in ('DELTA_PLUS_VAR', 'DELTA_MINUS_VAR'):
                if d not in thisvar.attrs:
                    continue
                deltaname = thisvar.attrs[d]
                if deltaname in self._varinfo:
                    continue
                self._varinfo[deltaname] \
                    = self._process_delta(thisname, deltaname)
                self._varinfo[deltaname]['vartype'] = 'D' #just like other deps
                self._varinfo[deltaname]['thisdim'] = dim
        for a in ('DELTA_PLUS_VAR', 'DELTA_MINUS_VAR'): #Process DELTA vars
            if not a in attrs:
                continue
            thisname = attrs[a]
            if thisname not in self._varinfo:
                #If DELTA_PLUS/DELTA_MINUS are same var, skip second one
                self._varinfo[thisname] \
                    = self._process_delta(self._name, thisname)
                self._varinfo[thisname]['vartype'] = 'U'

    def slice(self, dim, start=None, stop=None, step=None,
              single=False):
        """Slice on a single dimension

        Selects subset of a dimension to include in the output. Slicing
        is done with reference to the dimensions of the main variable and
        the corresponding dimensions of all other variables are sliced
        similarly. The first non-record dimension of the variable is always
        1; 0 is the record dimension (and is ignored for NRV variables).
        Multiple slices can be applied to select subsets of multiple
        dimensions; however, if one dimension is indexed multiple
        times, only the last one in the chain takes effect.

        Interpretation of the slice parameters is like normal Python slicing,
        including the ability to use negative values, etc.

        Passing in only a dimension "resets" the slice to include the
        entire dimension.

        Parameters
        ----------
        dim : int
            CDF dimension to slice on. This is the dimension as specified
            in the CDF (0-base for RV variables, 1-base for NRV) and does
            not change with successive slicing. Each dimension can only be
            sliced once.

        single : bool
            Treat ``start`` as a single index and return only that index
            (reducing dimensionality of the data by one.)

        start : int
            Index of first element of ``dim`` to include in the output.
            This can also be a sequence of indices to include, in which
            case ``stop`` and ``step`` must not be specified. This can be
            substantially slower than specifying ``stop`` and ``step``.

        stop : int
            Index of first element of ``dim`` to exclude from the output.

        step : int
            Increment between elements to include in the output.

        Returns
        -------
        VarBundle
            This bundle, for method chaining. This is not a copy: the
            original object is updated.

        Examples
        --------
        See the `VarBundle` examples for creating output from
        the slices.

        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
        >>> #Select index 2 from axis 1
        >>> b.slice(1, 2, single=True)
        >>> #Select from index 5 to end for axis 2, keeping index 2 from axis 1
        >>> b.slice(2, 5)
        >>> #Select 10 through 15 on axis 2, but all of axis 1
        >>> b.slice(1).slice(2, 10, 15)
        >>> #Select just record 5 and 10
        >>> b.slice(2).slice(0, [5, 10])
        >>> infile.close()
        """
        if single and (self._summed[dim] or self._mean[dim]):
            raise ValueError('Cannot sum/average on a single-element slice.')
        self._degenerate[dim] = single

        fancyidx = (stop is None and step is None and numpy.ndim(start) != 0)
        sl = slice(None, None, None) if fancyidx else slice(start, stop, step)
        for v in self._varinfo.values():
            if not dim in v['dims']:
                continue #This "main" var dimension isn't in this var
            idx = v['dims'].index(dim)
            #The slice to perform on read
            v['slice'][idx] = start if single else sl
            #And the slice to perform after the fact
            if fancyidx:
                v['postidx'][idx] = start
        return self

    def sum(self, dim):
        """Sum across a dimension.

        Total the main variable of the bundle across the given dimension.
        That dimension disappears from the output and dependencies
        (including their uncertainties) are assumed to be constant across
        the summed dimension. The uncertainty of the main variable, if
        any, is appropriately propagated (quadrature sum.)

        An invalid value for any element summed over will result in a fill
        value on the output. This does not work well for variables that
        define multiple VALIDMIN/VALIDMAX based on position within a
        dimension; the smallest VALIDMIN/largest VALIDMAX rather than the
        position-specific value.

        Summing occurs after slicing, to allow summing of a subset of
        a dimension. A single element slice (which removes the dimension)
        is incompatible with summing over that dimension.

        There is not currently a way to "undo" a sum; create a new
        bundle instead.

        Parameters
        ----------
        dim : int
            CDF dimension to total. This is the dimension as specified
            in the CDF (0-base for RV variables, 1-base for NRV) and does
            not change with successive slicing or summing. This must be a
            positive number (no support for e.g. -1 for last dimension.)

        Returns
        -------
        VarBundle
            This bundle, for method chaining. This is not a copy: the
            original object is updated.

        Examples
        --------
        See the `VarBundle` examples for creating output.

        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['Counts_P'])
        >>> #Total over dimension 1 (pitch angle)
        >>> b.sum(1)
        >>> #Get a new bundle (without the previous sum)
        >>> b = spacepy.pycdf.istp.VarBundle(infile['Counts_P'])
        >>> #Total over first 10 elements of dimension 2 (energy bins)
        >>> b.slice(2, 0, 10).sum(2)
        >>> infile.close()
        """
        if self._degenerate[dim]:
            raise ValueError('Cannot sum on a single-element slice.')
        if self._mean[dim]:
            raise ValueError('Cannot sum and take mean of same dimension.')
        self._summed[dim] = True
        return self

    def mean(self, dim):
        """Take the mean of a dimension.

        Take mean of the main variable of the bundle across the given
        dimension. That dimension disappears from the output and dependencies
        (including their uncertainties) are assumed to be constant across
        the summed dimension. The uncertainty of the main variable, if
        any, is appropriately propagated.

        Invalid values are excluded fromthe mean. This does not work well
        for variables that define multiple VALIDMIN/VALIDMAX based on
        position within a dimension; the smallest VALIDMIN/largest VALIDMAX
        rather than the position-specific value.

        Averaging occurs after slicing, to allow averaging of a subset of
        a dimension. A single element slice (which removes the dimension)
        is incompatible with averaging over that dimension.

        There is not currently a way to "undo" a mean; create a new
        bundle instead.

        Parameters
        ----------
        dim : int
            CDF dimension to average. This is the dimension as specified
            in the CDF (0-base for RV variables, 1-base for NRV) and does
            not change with successive slicing or summing. This must be a
            positive number (no support for e.g. -1 for last dimension.)

        Returns
        -------
        VarBundle
            This bundle, for method chaining. This is not a copy: the
            original object is updated.

        Examples
        --------
        See the `VarBundle` examples for creating output.

        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['Counts_P'])
        >>> #Average over dimension 1 (pitch angle)
        >>> b.mean(1)
        >>> #Get a new bundle (without the previous sum)
        >>> b = spacepy.pycdf.istp.VarBundle(infile['Counts_P'])
        >>> #Average over first 10 elements of dimension 2 (energy bins)
        >>> b.slice(2, 0, 10).mean(2)
        >>> infile.close()
        """
        if self._degenerate[dim]:
            raise ValueError('Cannot average on a single-element slice.')
        if self._summed[dim]:
            raise ValueError('Cannot sum and take mean of same dimension.')
        self._mean[dim] = True
        return self

    def _tokeep(self):
        """Determine which variables to keep for output

        Dependencies for dimensions which disappear after slicing, and
        other variables that they depend on, shouldn't be included in
        the output

        Returns
        -------
        list of str
            Names of variables to include in the output.
        """
        #What dims of main var disappear?
        deleted = [i for i in range(len(self._degenerate))
                   if any((self._degenerate[i], self._summed[i],
                           self._mean[i]))]
        return [v for v, i in self._varinfo.items()
                if i.get('thisdim', None) not in deleted]

    def _same(self, newvar, invar, rv, dv, dims, data):
        """Checks if an existing variable matches a proposed new variable

        Does not compare DEPEND and LABL_PTR attributes (those are handled
        later.)

        Parameters
        ----------
        newvar : `~.pycdf.Var`
            Existing variable to compare to requirements

        invar : : class:`~.pycdf.Var`
            Variable to use as reference for attributes, RV, CDF type,
            number of elements.

        rv : bool
            Is the new variable record-varying

        dv : list of bool
            Data variance for each dimension.

        dims : list of int
            Size of each dimension.

        data : `~numpy.ndarray`
            Data that should be in the variable.

        Returns
        -------
        bool
            True if the existing variable is the same; False if not.
        """
        # CDF output only checks
        if hasattr(newvar, 'type'):
            if newvar.rv() != rv or newvar.dv() != dv:
                return False
            if hasattr(invar, 'type') and newvar.type() != invar.type():
                return False
            if hasattr(invar, 'nelems') and newvar.nelems() != invar.nelems():
                return False
        #Check basic type, dimensions, etc.
        if newvar.dtype != invar.dtype\
            or len(dims) != (len(newvar.shape) - rv) \
            or list(dims) != list(newvar.shape[rv:]):
            return False
        ia = invar.attrs
        na = newvar.attrs
        for a in ia:
            if a.startswith(('DEPEND_', 'LABL_PTR_')) \
               or a == 'FIELDNAM':
                #depends/LABL PTR shift around, and FIELDNAM may change,
                #so test outside of this function.
                pass
            if not a in na or not numpy.array_equal(ia[a], na[a]):
                return False
            # CDF input *and* output only
            if hasattr(na, 'type') and hasattr(ia, 'type')\
               and ia.type(a) != na.type(a):
                return False
        #Finally check the data
        return (data == newvar[...]).all()

    def _namemap(self, suffix=None):
        """Map old variable names to new

        Helper for `output` that maps the variable name in the
        input CDF to variable name in the output CDF.

        Parameters
        ----------
        suffix : str
            String to append to name of variables that are changed
            from input to output.

        Returns
        -------
        dict
            Keyed by name in input, values are name in the output. No
            entry for names that don't change.
        """
        namemap = {}
        if suffix is not None:
            for vname, vinfo in self._varinfo.items():
                if vinfo['vartype'] in ('M', 'U'):
                    namemap[vname] = vname + suffix
                else: #Dependency. If any slice/sum, it's changed
                    if any([any((self._summed[d], self._mean[d]))
                            for d in vinfo['dims']]) \
                    or any([s != slice(None) for s in itertools.chain(
                        vinfo['slice'], vinfo['postidx'])]):
                        namemap[vname] = vname + suffix
        return namemap

    def _sum_avg(self, data, invar, vinfo, degen, summed, averaged):
        """Sum/average data

        Helper for `output` that performs summing and averaging
        of the data for a single variable. Note dimensionality of all
        input is before the removal of degenerate dimensions
        (this function does the translation using ``degen``), and it is
        by dimension not axis (so NRV variables have a vestigial 0th
        dimension that is not interpreted.)

        Parameters
        ----------
        data : `numpy.ndarray`
            Data as read from input CDF and properly sliced.

        invar : `~.pycdf.Var`
            CDF input variable from which ``data`` was read.

        vinfo : dict
            Value from instance variable ``_varinfo`` for this variable.

        degen : list of bool
            For each dimension of this variable, whether the dimension
            is degenerate (i.e. already gone at this point.)

        summed : list of bool
            For each dimension of this variable, whether the dimension
            should be summed over.

        averaged : list of bool
            For each dimension of this variable, whether the dimension
            should be averaged over.

        Returns
        -------
        `numpy.ndarray`
            Data summed/averaged over dimensions according to ``summed``
            and ``averaged`` inputs.
        """
        #Correction for NRV variables in the mapping between dim and axis
        nrv = int(not vinfo['rv'])
        #Degenerate slices have already been removed, so need
        #a map from old dim numbers to new ones. Note removing
        #the record dimension does not shift other dims!
        newdims = [None if degen[i] else i - sum(degen[1:i])
                   for i in range(len(degen))]
        #Axis numbers to sum, with degenerate removed
        #(NRV means dim 0 is axis 1, so correct for that,
        #and also don't do any actions on dim 0 for NRV)
        summe = [newdims[i] - nrv
                 for i in range(nrv, len(summed)) #old dim
                 if newdims[i] is not None and (summed[i] or averaged[i])]
        avgme = [newdims[i] - nrv
                 for i in range(nrv, len(averaged)) #old dim
                 if newdims[i] is not None and averaged[i]]
        #Sum over axes in reverse order so axis renumbering
        #doesn't affect future sums
        a = invar.attrs
        for ax in summe[::-1]:
            if vinfo['vartype'] == 'D':
                #If sum over DEPEND, must be constant over axis
                data = data.take(0, axis=ax)
                continue
            invalid = numpy.isclose(data, a['FILLVAL'])
            if 'VALIDMIN' in a:
                invalid = numpy.logical_or(
                    invalid, data < numpy.min(a['VALIDMIN']))
            if 'VALIDMAX' in a:
                invalid = numpy.logical_or(
                    invalid, data > numpy.max(a['VALIDMAX']))
            data[invalid] = 0 #avoids warning and helps with mean
            if vinfo['vartype'] == 'M':
                data = data.sum(axis=ax)
            elif vinfo['vartype'] == 'U': #propagate error
                data = numpy.sqrt((data ** 2).sum(axis=ax))
            else: #Should not happen
                raise ValueError('Bad summation type.')
            if ax in avgme: #divide out
                count = numpy.sum(~invalid, axis=ax, dtype=data.dtype)
                invalid = (count == 0)
                count[invalid] = 1 #avoid warning
                data = data / count
            else: #Sum, so any fill on axis means value is fill
                invalid = invalid.max(axis=ax)
            data[invalid] = a['FILLVAL']
        return data

    def _repoint_depend(self, invar, newvar, preexist, namemap, degen):
        """Change DEPEND for new dimensionality of one variable.

        Slicing/summing might change variable dimensionality and thus
        the relationship with its DEPENDs, and the DEPENDs themselves
        may have a new name. This updates the DEPEND attributes for
        these changes, or verifies they are correct if the output
        variable already exists.

        Parameters
        ----------
        invar : `~.pycdf.Var`
            The input variable (opened in raw mode).

        newvar : `~.pycdf.Var`
            The output variable (opened in raw mode).

        preexist : bool
            True if ``newvar`` existed and doing a consistency check;
            False if ``newvar`` was newly created and should be edited.

        namemap : dict
            Map from name in input variable (key) to name in output
            variable (value). No entry if name didn't change.

        degen : list of bool
            For each dimension of this variable, whether the dimension
            is degenerate (i.e. already gone at this point.) This
            includes any degeneracy from summing/averaging as well as
            slicing.


        """
        #Index by old dim; returns the new dim (None if went away)
        #Note slicing away DEPEND_0 (record dimension) does NOT change
        #subsequent depends!
        newdims = [None if degen[i] else i - sum(degen[1:i])
                   for i in range(len(degen))]
        for a in list(newvar.attrs.keys()): #Editing in loop!
            #Handle a suffixed DELTA if necessary
            if a.startswith('DELTA_'):
                olddelta = invar.attrs[a]
                if isinstance(olddelta, bytes):
                    olddelta = olddelta.decode('ascii')
                newvar.attrs[a] = namemap.get(olddelta, olddelta)
                continue
            if not a.startswith(('DEPEND_', 'LABL_PTR_')):
                continue
            newdim = int(a.split('_')[-1])
            oldval = None #Sentinel value
            if newdim in newdims: #An old value that belongs in this dim
                olddim = newdims.index(newdim)
                old_a = '{}_{}'.format('_'.join(a.split('_')[:-1]),
                                       olddim)
                oldval = invar.attrs.get(old_a, None)
                if isinstance(oldval, bytes):
                    oldval = oldval.decode('ascii')
            if oldval is not None:
                #Check for variable renaming from the input to output
                newval = namemap.get(oldval, oldval)
                if preexist:
                    existingval = newvar.attrs[a]
                    if isinstance(existingval, bytes):
                        existingval = existingval.decode('ascii')
                    if existingval != newval:
                        raise RuntimeError(
                            'Incompatible {} already exists in output.'
                            .format(newvar.name()))
                else:
                    newvar.attrs[a] = newval
            else:
                #Either there's no corresponding old dim, or it didn't have
                #a DEPEND. Either way, shouldn't be a DEPEND in the new dim.
                if preexist:
                    if a in newvar.attrs:
                        raise RuntimeError(
                            'Incompatible {} already exists in output.'
                            .format(newvar.name()))
                else:
                    del newvar.attrs[a]

    def _outshape(self, vname):
        """Calculate shape of the variable on output

        Parameters
        ----------
        vname : str
            Name of the variable to check the shape of.

        Returns
        -------
        tuple
            The shape of the variable after all slicing, etc. applied, or
            None of the variable is not included in output.
        """
        if vname not in self._tokeep():
            return None
        vinfo = self._varinfo[vname]
        invar = self.cdf[vname]
        rv = vinfo['rv']
        shape = invar.shape
        sl = vinfo['slice']
        postidx = vinfo['postidx']
        #no dimension has BOTH a slice and a postindex, so combine
        slices = [pi if s == slice(None, None, None) else s
                  for s, pi in zip(sl, postidx)]
        #And any dim that is summed/averaged is degenerate, so
        #slice with a single index to make it go away
        for d in vinfo['dims']:
            if self._summed[d] or self._mean[d]:
                slices[d] = 0
        if not rv: #Remove record dimension
            slices = slices[1:]
        #Make a fake array the size of the input, and slice it
        return numpy.empty(shape=shape)[tuple(slices)].shape

    def variables(self):
        """Description of variable output from the bundle

        Provides information describing the variables output
        from the bundle

        Returns
        -------
        list
            Each element is a list-of-tuples. The list corresponds to a
            dimension of the master var: first the master var itself, then the
            uncertainties and labels associated with each dimension. Each
            element of these sublists is then a tuple of variable name and
            shape on the output (itself a tuple). If a variable isn't
            included in the output (sliced away), its shape will be ``None``.

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
        >>> b.slice(1, 2, single=True).variables()
        [[('FPDU', (100, 72))],
         [('Epoch_Ion', (100,)), ('Epoch_Ion_DELTA', (100,))],
         [('PITCH_ANGLE', None), ('Pitch_LABL', None)],
         [('HOPE_ENERGY_Ion', (100, 72)),
          ('ENERGY_Ion_DELTA', (100, 72)),
          ('Energy_LABL', (72,))]]
        """
        #List of every variable in each dimension
        v_by_dim = functools.reduce(
            lambda x, vname:
            x[self._varinfo[vname].get('thisdim', None)].append(vname) or x,
            self._varinfo.keys(), collections.defaultdict(list))
        for l in v_by_dim.values():
            l.sort(key=lambda x: (self._varinfo[x]['sortorder'], x))
        variables = [[(v, self._outshape(v))
                      for v in v_by_dim.get(None, [])]]
        vi = self._varinfo[self._name]
        for dim in vi['dims']:
            variables.append([
                (v, self._outshape(v)) for v in v_by_dim.get(dim, [])])
        return variables

    def operations(self):
        """Operations of this bundle

        Provides information describing the operations this bundle
        would perform.

        Returns
        -------
        list
            Each element is a tuple: first element is a string with
            the name of the operation (i.e. method of
            `VarBundle`), next is also a tuple of positional
            arguments, and finally a dict of keyword arguments.

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
        >>> b.slice(1, 2, single=True).operations()
          [('slice', (1, 2), {'single': True})]
        >>> #Apply same operations to a different variable
        >>> b2 = spacepy.pycdf.istp.VarBundle(infile['FEDU'])
        >>> for op, args, kwargs in b2.operations():
        ...     getattr(b2, op)(*args, **kwargs)
        """
        ops = []
        vi = self._varinfo[self._name]
        for dim in vi['dims']:
            sl = vi['slice'][dim]
            postidx = vi['postidx'][dim]
            if sl != slice(None, None, None): #simple slice
                if isinstance(sl, slice): #slice
                    ops.append((
                        'slice',
                        tuple((s for s in (dim, sl.start, sl.stop, sl.step)
                               if s is not None)),
                        {}))
                else: #single index
                    ops.append(('slice', (dim, sl), {'single': True}))
            elif postidx != slice(None, None, None): #fancy index
                ops.append(('slice', (dim, postidx,), {}))
            for v, name in zip((self._mean, self._summed), ('mean', 'sum')):
                if v[dim]:
                    ops.append((name, (dim,), {}))
        return ops

    def output(self, output, suffix=None):
        """Output the variables as modified

        Parameters
        ----------
        output : `~.pycdf.CDF`,  `~.datamodel.SpaceData`
            Output container to receive the new data, may be an open CDF
            file or a SpaceData.

        suffix : str
            Suffix to append to the name of any variables that are changed
            for the output. This allows the output to contain multiple
            variables derived from the same input variable. The main variable
            and its DELTA variables will always have the suffix applied.
            Any dependencies will have the suffix applied only if they have
            changed from the input CDF (e.g. from slicing.)

        Returns
        -------
        VarBundle
            This bundle, for method chaining.

        See Also
        --------
        toSpaceData

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
        >>> outfile = spacepy.pycdf.CDF('output.cdf', create=True)
        >>> #Output the low energy half in one variable
        >>> b.slice(2, 0, 36).output(outfile, suffix='_LoE')
        >>> #And the high energy half in another variable
        >>> b.slice(2, 36, 72).output(outfile, suffix='_HiE')
        >>> outfile.close()
        >>> infile.close()
        """
        tokeep = self._tokeep()
        namemap = self._namemap(suffix)
        for vname in tokeep:
            vinfo = self._varinfo[vname]
            #Dim of main var that depends on this (None if main var or delta)
            maindim = vinfo.get('thisdim', None)
            #Degeneracy of dimensions in this variable's "frame"
            degen = [self._degenerate[d] for d in vinfo['dims']]
            #And whether the dim was summed
            summed = [self._summed[d] for d in vinfo['dims']]
            #And averaged
            averaged = [self._mean[d] for d in vinfo['dims']]
            # Raw data for CDF input *and* output only
            invar = self.cdf.raw_var(vname) if hasattr(output, 'raw_var')\
                    and hasattr(self.cdf, 'raw_var') else self.cdf[vname]
            sl = vinfo['slice'] #including 0th dim
            postidx = vinfo['postidx']
            #Dimension size/variance for original variable
            #(0 index is CDF dimension 1)
            dv = self._varinfo[vname]['dv']
            rv = self._varinfo[vname]['rv']  #and record variance
            #Scrub degenerate dimensions from the post-indexing
            #(record is never degenerate)
            postidx = [postidx[i] for i in range(len(postidx))
                       if not degen[i]]

            #Now get the data, and sum/average it
            if not rv: #Remove fake record dimension
                sl = sl[1:]
                postidx = postidx[1:]
            #Forces array scalars, makes the rest work better
            data = numpy.asanyarray(invar.__getitem__(tuple(sl)))
            if postidx:
                data = data[tuple(postidx)]
            data = self._sum_avg(data, invar, vinfo, degen, summed, averaged)
            #Summed/averaged dimensions are now also degenerate
            degen = [max(v) for v in zip(degen, summed, averaged)]

            #Get shape of output variable from actual data
            dims = data.shape
            #Raw Epoch16 have a trailing (2,)
            if hasattr(invar, 'type')\
               and invar.type() == spacepy.pycdf.const.CDF_EPOCH16.value:
                dims = dims[:-1]
            #Cut out any degenerate dimensions from DV (skipping record dim)
            dv = [dv[i] for i in range(len(dv)) if not degen[i + 1]]
            #Change record variance for the output if sliced away 0th
            if rv and degen[0]:
                rv = False
            if rv: #remove record dimension from size IF output is RV
                dims = dims[1:]

            #Rename the variable if necessary
            outname = namemap.get(vname, vname)
            if outname in output:
                preexist = True
                newvar = output.raw_var(outname) if hasattr(output, 'raw_var')\
                         and hasattr(self.cdf, 'raw_var') else output[outname]
                if not self._same(newvar, invar, rv, dv, dims, data):
                    raise RuntimeError(
                        'Incompatible {} already exists in output.'
                        .format(outname))
            else:
                preexist = False
                if hasattr(output, 'new'):
                    t = invar.type() if hasattr(invar, 'type') else None
                    try:
                        compress, compress_param = invar.compress()
                    except (TypeError, AttributeError):
                        # arrays have a different "compress"
                        compress, compress_param = None, None
                    ne = invar.nelems() if hasattr(invar, 'nelems') else None
                    newvar = output.new(
                        outname, data=data,
                        type=t, recVary=rv,
                        dimVarys=dv, dims=dims,
                        n_elements=ne,
                        compress=compress, compress_param=compress_param)
                    newvar.attrs.clone(invar.attrs)
                else:
                    newvar = spacepy.dmarray(data, attrs=invar.attrs.copy())
                    output[outname] = newvar
                if vname != outname: #renamed
                    newvar.attrs['FIELDNAM'] = outname

            self._repoint_depend(invar, newvar, preexist, namemap, degen)
        return self

    def toSpaceData(self, suffix=None):
        """Return variables, as modified.

        Convenience function to call `output` on a new
        `~.datamodel.SpaceData` and return it.

        Parameters
        ----------
        suffix : str
            Appended to the name of variables changed on output; see
            `output` for details.

        Returns
        -------
        `.datamodel.SpaceData`
            Data read from input and processed according to the defined
            operations.

        See Also
        --------
        output

        Examples
        --------
        >>> import spacepy.pycdf
        >>> import spacepy.pycdf.istp
        >>> infile = spacepy.pycdf.CDF('rbspa_rel04_ect-hope-PA-L3_20121201_v7.1.0.cdf')
        >>> b = spacepy.pycdf.istp.VarBundle(infile['FPDU'])
        >>> data = b.slice(1, 2, single=True).toSpaceData()
        >>> infile.close()
        >>> data.tree()
        +
        |____ENERGY_Ion_DELTA
        |____Energy_LABL
        |____Epoch_Ion
        |____Epoch_Ion_DELTA
        |____FPDU
        |____HOPE_ENERGY_Ion
        """
        sd = spacepy.SpaceData()
        self.output(sd, suffix=suffix)
        return sd

    @staticmethod
    def _vtype(v):
        """String representation of type of a variable

        Parameters
        ----------
        v
            Open CDF variable, numpy array, or similar

        Returns
        -------
        str
            String representation of type of ``v``, either as CDF type
            or numpy type
        """
        # Kludge, but assumes main CDF code gets it right
        res = str(v).split(' ')[0]
        if res.startswith('CDF_'):
            return res
        return str(v.dtype)

    def __str__(self):
        """String representation of the bundle

        Returns a string representation of the bundle, which is all the
        variables that are involved on the input. Variables which are
        not included on the output are in []

        Returns
        -------
        str
            Brief string description of the bundle.
        """
        return '\n'.join([
            '{}{}: {} {}{}'.format(
                ' ' * 4 if self._varinfo[vname]['sortorder'] > 1 else '',
                vname,
                self._vtype(self.cdf[vname]),
                str(list(shape)) if shape is not None else '---',
                #RV vars always have dim 0 as axis 0, so they become
                #NRV iff dim 0 of the main var goes away
                ' NRV' if shape is not None
                and (not self._varinfo[vname]['rv'] or max(
                    self._degenerate[0], self._summed[0], self._mean[0]))
                else ''
            )
            for dimvars in self.variables() for vname, shape in dimvars])

    def __repr__(self):
        """Representation of bundle

        Cannot return anything that can be evaluated to create a copy
        of the CDF, so this is just the informal str representation in
        angle brackets.

        Returns
        -------
        str
            Informal representation of bundle contents.
        """
        return '<VarBundle:\n{}\n>'.format(str(self))
