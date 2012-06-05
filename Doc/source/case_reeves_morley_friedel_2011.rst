*************************************************
Paulikas and Blake revisited (Reeves et al. 2011)
*************************************************

This case study reproduces the figures of Reeves et al. (2011),
"On the relationship between relativistic electron flux and solar wind
velocity: Paulikas and Blake revisisted"
(`doi:10.1029/2010JA015735 <http://dx.doi.org/10.1029/2010JA015735>`_).

Setup
=====
Create a directory to hold files for this case study. Within this
directory, create subdirectories ``code``, ``data``, and
``plots``. (Using version control on the code directory is
recommended; the SpacePy team uses `git
<http://git-scm.com/documentation>`_.)

Obtaining energetic particle data
=================================
We require the 1.8-3.5 MeV electron flux from the LANL-GEO ESP
detector, available in the paper's `auxiliary material
<http://www.agu.org/journals/ja/ja1102/2010JA015735/supplement.shtml>`_. The
ESP data are in Data Set S1. Save this file to the ``data`` directory;
the filename is assumed to be ``2010ja015735-ds01.txt``.

The data file was corrupted on upload to AGU, and the code to fix it
is non-trivial, so this is a good chance to learn how to run someone
else's code. (:ref:`appendix` has step-by-step information on each
portion of this process.) Copy all of the following and paste it into
a file called ``fix_esp_dta.py`` in the ``code`` directory.

.. code-block:: python

    import os.path


    datadir = os.path.join('..', 'data')
    in_name = os.path.join(datadir, '2010ja015735-ds01.txt')
    out_name = os.path.join(datadir, '2010ja015735-ds01_FIXED.txt')
    infile = open(in_name, 'r')
    outfile = open(out_name, 'w')
    data = infile.read()
    infile.close()

    data = data.replace('\r', '\n')
    data = data.replace('\n\n', '\n')
    data = data.split('\n')

    for i in range(15):
        outfile.write(data.pop(0) + '\n')
    oldline = None
    for line in data:
        if line[0:2] in ['19', '20', '2']:
            if not oldline is None:
                outfile.write(oldline + '\n')
            oldline = line
        else:
            oldline += line
    outfile.write(oldline + '\n')
    outfile.close()

Now this script can be run with ``python fix_esp_data.py``. It should
create a file called ``2010ja015735-ds01_FIXED.txt`` in the ``data``
directory.

Edit this with your favorite text editor and fix the fluxes for 1993-4-13,
1996-4-24, 2000-1-7, 2009-9-10. Each of these is missing the decimal point
in the flux.

Data fixed, we can load and begin examining the data.  Change to the
``code`` directory and start your Python interpreter. (`IPython
<http://ipython.org/>`_ is recommended, but not required.)

In the following examples, do not type the leading ``>>>``; this is
the Python interpreter prompt. IPython has a different prompt that
looks like ``In [1]``.

>>> import os.path
>>> datadir = os.path.join('..', 'data')
>>> print(datadir)
../data

The first line imports the :py:mod:`os.path` module from the Python
standard library. Python has a huge `standard library
<http://docs.python.org/library/index.html>`_. To keep this code
organized, it is divided into many modules, and a module must be
imported before it can be used. (The `Python module of the week
<http://www.doughellmann.com/PyMOTW/>`_ is a great way to explore the
standard library.)

The second line makes a variable, ``datadir``, which will contain the
path of the data directory. The :py:func:`os.path.join` function
provides a portable way of "gluing" together directories in a path,
and will use backslashes on Windows and forward slashes on Unix. The
third line then prints out the value of this variable for
confirmation; note this is a Unix system.

Note that string constants in Python can use single or double quotes;
we could just as well have written:

>>> datadir = os.path.join("..", "data")

or even:

>>> datadir = os.path.join('..', "data")

The full path can also be used (and this is a better case for using a
variable.) For example, I am preparing this example in a directory
``reeves_morley_friedel_2011`` in my home directory, so I could use:

>>> datadir = os.path.join('home', 'jniehof', 'reeves_morley_friedel_2011',
...                        'data')

This very long line can be typed across two lines in Python, and
because the line break happens within parentheses, a line continuation
character is not required.

Returning to reading the ESP data file:

>>> fname = os.path.join(datadir, '2010ja015735-ds01_FIXED.txt')

creates a variable holding the full path to the fixed file.

>>> import numpy

The import statement imports any installed `module <http://docs.python.org/tutorial/modules.html>`_, just as if it were in the standard library. Here we import the very useful :mod:`numpy` module, which is a prerequisite for SpacePy and useful in its own right.

>>> esp_fluxes = numpy.loadtxt(fname, skiprows=15, usecols=[2])

:func:`~numpy.loadtxt` makes it easy to load data from a file into a numpy :class:`~numpy.ndarray`, a very useful data container. ``skiprows`` skips the header information, and specifying only column 2 with ``usecols`` will only load the fluxes for 1.8-3.5MeV. We only load the fluxes at this point because they can be represented as floats, which numpy arrays store very efficiently.

>>> import datetime

The :mod:`datetime` module provides Python objects which can manipulate dates and times and have some understanding of the meanings of dates, making for easy comparisons between dates, date arithmetic, and other useful features.

>>> convert = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d')

This line sets up a converter to be used later. :meth:`~datetime.datetime.strptime` creates a :class:`~datetime.datetime` from a string, given a format definition (here specified as year-month-day). So:

>>> print(datetime.datetime.strptime('2010-01-02', '%Y-%m-%d'))
2010-01-02 00:00:00

`lambda
<http://docs.python.org/tutorial/controlflow.html#lambda-forms>`_ is a
simple shortcut for a one-liner function; wherever ``convert(x)`` is
used after the definition, it functions like
``datetime.datetime.strptime(x, '%Y-%m-%d')``. This makes it easier to
parse a date string without specifying the format all the time:

>>> print(convert('2010-01-02'))

This converter can be used with :func:`~numpy.loadtxt`:

>>> esp_times = numpy.loadtxt(fname, skiprows=15, usecols=[0,],
...                           converters={0: convert}, dtype=numpy.object)

The ``converters`` option takes a Python `dictionary
<http://docs.python.org/tutorial/datastructures.html#dictionaries>`_. The
default `dtype
<http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html>`_ is
float, which cannot store datetimes; using ``numpy.object``
allows storage of any Python object.

Since it would be useful to be able to load the data without typing so
many lines, create a file called ``common.py`` in the ``code``
directory with the following contents:

.. code-block:: python

    import datetime
    import os.path

    import numpy


    datadir = os.path.join('..', 'data')

    def load_esp():
        fname = os.path.join(datadir, '2010ja015735-ds01_FIXED.txt')
        esp_fluxes = numpy.loadtxt(fname, skiprows=15, usecols=[2])
        convert = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d')
        esp_times = numpy.loadtxt(fname, skiprows=15, usecols=[0,],
                                  converters={0: convert}, dtype=numpy.object)
        return (esp_times, esp_fluxes)

All needed imports are at the top of the file, with one blank line
between standard library imports and other imports and two blank lines
after them. ``datadir`` is defined as a global variable, outside of
the function (but notice that it is available to the ``load_esp``
function.)

The rest of the file defines a `function
<http://docs.python.org/tutorial/controlflow.html#defining-functions>`_
which returns the dates and fluxes in a `tuple
<http://docs.python.org/tutorial/datastructures.html#tuples-and-sequences>`_. The
next section shows how to use this function.

Solar Wind data and averaging
=============================

The top panel of figure 1 shows the ESP fluxes overplotted with the
solar wind velocity. Fortunately, the :mod:`~spacepy.omni` module of
SpacePy provides an interface to the hourly solar wind dataset,
OMNI. The data are stored in a `dictionary
<http://docs.python.org/tutorial/datastructures.html#dictionaries>`_
called ``omnidata``, which we will access directly since we do not
need the interpolation functions of :func:`~spacepy.omni.get_omni`:

>>> import spacepy.omni
>>> vsw = spacepy.omni.omnidata['velo']
>>> vsw_times = spacepy.omni.omnidata['UTC']

We'll also load the esp data:

>>> import common
>>> esp_times, esp_flux = common.load_esp()

Even though we have not installed ``common.py``, the ``import``
statement finds it because it is in the current directory.

``load_esp`` returns a `tuple
<http://docs.python.org/release/2.6.7/tutorial/datastructures.html#tuples-and-sequences>`_,
which can be *unpacked* into separate variables.

Now we need to produce 27-day running averages of both the flux and
the solar wind speed. Fortunately there are no gaps in the time
series:

>>> import numpy
>>> d = numpy.diff(vsw_times)
>>> print(d.min())
1:00:00
>>> print(d.max())
1:00:00
>>> d = numpy.diff(esp_times)
>>> print(d.min())
1 day, 0:00:00
>>> print(d.max())
1 day, 0:00:00

:func:`numpy.diff` returns the difference between every element of an
array and the previous element. :meth:`~numpy.ndarray.min` and
:meth:`~numpy.ndarray.max` do exactly what they sound like. So this
code confirms that every time in the vsw data is on a continuous one
hour cadence, and the ESP data is on a continuous one day cadence.

>>> import scipy
>>> esp_flux_av = numpy.empty(shape=esp_flux.shape, dtype=esp_flux.dtype)
>>> for i in range(len(esp_flux_av)):
...     esp_flux_av[i] = scipy.stats.nanmean(esp_flux[i - 13:i + 14])

:func:`numpy.empty` creates an empty array, taking the ``shape`` and
``dtype`` from the ``esp_flux`` array. ``empty`` does not initalize
the data in the array, so it is essentially random junk; use
:func:`~numpy.zeros` to create an array filled with zeros.

:func:`len` returns the length of an array, and :func:`range` then
iterates over each number from 0 to length minus 1, i.e. the entire
array. Each element is then set to a 27-day average: from 13 days
before a day's measurement through 13 days after. (Python slices do
not include the last element listed; they are half-open). Note that
these slices can happily run off the end or the beginning of the
``esp_flux`` array.

:func:`~scipy.stats.stats.nanmean` takes the mean of a numpy array,
but skips any elements with a value of "not a number" (nan), which is
often used for fill.  (This is our first exposure to the :mod:`scipy`
module.)

The solar wind data covers from 1963, whereas the ESP data starts in
1989.  Although for proper averaging we want to keep some solar wind
data "off the end" of the ESP data, 35 years is a bit much. So let's
cut out the solar wind data from before 1989:

>>> import bisect
>>> import datetime
>>> idx = bisect.bisect_left(vsw_times, datetime.datetime(1989, 1, 1))
>>> vsw_times = vsw_times[idx:]
>>> vsw = vsw[idx:]

:mod:`bisect` provides fast functions for searching in sorted data;
:func:`~bisect.bisect_left` is roughly a find-the-position-of function.
Having found the position of the start of 1989, we then keep times
from then on (specifying a start index without a stop index in Python
means "from start to end of the list.") Note that, although ``bisect``
is meant to work on lists, it works fine on numpy arrays; this is a
common feature of Python known as
`duck typing <http://en.wikipedia.org/wiki/Duck_typing#In_Python>`_.

For the solar wind averaging, the times need to cover the 24 * 13.5 = 324
hours previous, and 324 hours following (non-inclusive). There is also a 
more efficient way than using an explicit loop:

>>> vsw_av = numpy.fromiter((scipy.stats.nanmean(vsw[i - 324:i + 324])
...                         for i in range(len(vsw))),
...                         count=len(vsw), dtype=vsw.dtype)

:func:`~numpy.fromiter` makes a numpy array from an `iterator
<http://docs.python.org/library/stdtypes.html#iterator-types>`_, which
is like a list except that it holds information on generating each
element in a sequence rather than creating the entire
sequence. ``count`` provides numpy with the number of elements in the
output (so it can make the entire array at once); ``dtype`` here is
just copied from the input.

The type of iterator used here is a `generator expression
<http://www.python.org/dev/peps/pep-0289/>`_, closely related to a
`list comprehension
<http://docs.python.org/tutorial/datastructures.html#list-comprehensions>`_.
These are among the most powerful and most difficult to understand
concepts in Python. An illustrative, although not useful, example:

>>> for i in (x + 1 for x in range(10)):
...     print(i)

Here ``(x + 1 for x in range(10))`` is a generator expression that
creates an iterator, which will return the numbers 1 through 10. At no
point is the complete list of all numbers constructed, saving memory.

In our calculation of ``esp_flux_av``, we created an explicit loop in
Python. The generator expression used to compute ``vsw_av`` has no
explicit loop, and the actual looping is handled in (much faster)
compiled C code.

.. _appendix:

Appendix: Fixing the ESP data file
=============================
This appendix provides a detailed explanation of the script that fixes
the ESP data file.

First set up a variable to hold the location of the data, as above:

>>> import os.path
>>> datadir = os.path.join('..', 'data')

Examining the data file, it is clear that something is odd: lines
appear to have been broken inappropriately; for example, the data for
1989-10-12 are split across two lines. So the first task is to fix
this file, first opening the original (broken) file and an output
(fixed) file:

>>> in_name = os.path.join(datadir, '2010ja015735-ds01.txt')
>>> out_name = os.path.join(datadir, '2010ja015735-ds01_FIXED.txt')
>>> infile = open(in_name, 'r')
>>> outfile = open(out_name, 'w')

These lines :func:`open` the original file for reading (``r``), and a
new file for writing (``w``). Note that opening a file for writing
will destroy any existing contents.

The file happens to contain a mixture of carriage returns and proper newlines, so to begin all the carriage returns need to be rewritten as newlines:

>>> data = infile.read()
>>> infile.close()
>>> data = data.replace('\r', '\n')
>>> data = data.replace('\n\n', '\n')

:meth:`~file.read` reads *all* data from the file at once, so this is
not recommended for large files. In this case it makes things
easier. Once the data are read, :meth:`~file.close` the file. Calling
the :meth:`~str.replace` method on ``data`` replaces all instances of
the first parameter (``'\r'``) with the second (``'\n'``). ``\r`` is
the special code indicating a carriage return; ``\n``, a newline. For
a literal backslash, use ``\\``. Once the carriage returns have been
replaced with newlines, a second round of replacement eliminates
duplicates.

Now that the line endings have been cleaned up, it's time to rejoin the erroneously split lines. First copy over the 15 lines of header verbatim:

>>> data = data.split('\n')
>>> for i in range(15):
...     outfile.write(data.pop(0) + '\n')

:meth:`~str.split` splits a string into a `list
<http://docs.python.org/tutorial/introduction.html#lists>`_, with the
split between elements happening wherever the provided parameter
occurs. A simple example:

>>> foo = 'a.b.c'.split('.')
>>> print(foo)
['a', 'b', 'c']

The splitting character is not present in the output.

The advantage of a list is that it makes it easy to access individual elements:
>>> print(foo[1])
b

The first element of a Python list is numbered zero.

:func:`range` returns a list of numbers, starting from 0, with the parameter specifying how many elements are in the list:

>>> print(range(5))
[0, 1, 2, 3, 4]

The last number is 4 (not 5 as might be expected), but there are 5
elements in the list.

The `for <http://docs.python.org/tutorial/controlflow.html#for-statements>`_ executes the following indented statement once for every element in the ``in`` list:

>>> for i in ['a', 'b', 'c']:
...     print i
a
b
c

Indentation is significant in Python! Normally indents are four spaces and the tab key will do the job. (In the above example, you may need to hit enter twice after the print statement, the second to terminate the indentation.)

`pop <http://docs.python.org/tutorial/datastructures.html#more-on-lists>`_ returns one element from a list, and deletes it from the list. Using ``0`` pops off the first element, and :meth:`~file.write` writes a string to a file. ``+`` can be used to concatenate two strings togther. Since :meth:`~str.split` removed the newlines, they need to be readded.

So this little block of code splits the data into a list on newlines and, repeating fifteen times, takes the first element of that list and writes it, with a newline, to the output. Now ``data`` contains only the actual lines of data.

>>> oldline = None
>>> for line in data:
...     if line[0:2] in ['19', '20', '2']:
...         if not oldline is None:
...             outfile.write(oldline + '\n')
...         oldline = line
...     else:
...         oldline += line
>>> outfile.write(oldline + '\n')
>>> outfile.close()

``None`` is a special Python value specifically indicating nothing;
it's used here to mark the first time around the loop.

``line[0:2]`` gets the first two characters in the string `line`, and
the ``in`` operator compares the resulting string to see if it is
present in the following list. This will return ``True`` if the line
begins with ``19`` or ``20``. The `if
<ttp://docs.python.org/tutorial/controlflow.html#if-statements>`_
statement executes the following indented block if the condition is
True. So, if this is True, the previous line probably ended properly
and it can be written out. First there is an additional check that
this isn't the first time around the loop, and then the *previous*
line (which we know ended cleanly) is written out. The currently-read
line then becomes the new "previous" line.

The ``2`` is a special case: if the line is less than two characters
long, ``line[0:2]`` will return the entire line, and it so happens
that these cases always correspond to the previous line being whole.

If this test fails, everything under ``else`` is executed. Here the
assumption is that the previous line didn't end cleanly and the
current line is actually a continuation of it, so the current line is
appended to the previous. ``a += b`` is a shortcut for ``a = a + b``.

Once the loop terminates, the last line is written out, and the file closed.
