*************************************************
Paulikas and Blake revisited (Reeves et al. 2011)
*************************************************

This case study reproduces the figures of Reeves et al. (2011),
"On the relationship between relativistic electron flux and solar wind
velocity: Paulikas and Blake revisited"
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
detector, available in the paper's auxiliary material (scroll down to "Supporting information" on the `paper's page <http://dx.doi.org/10.1029/2010JA015735>`_. The
ESP data are in Data Set S1. Save this file to the ``data`` directory;
the filename is assumed to be ``jgra20797-sup-0003-ds01.txt``.

The data file was corrupted on upload to AGU, and the code to fix it
is non-trivial, so this is a good chance to learn how to run someone
else's code. (:ref:`appendix` has step-by-step information on each
portion of this process.) Copy all of the following and paste it into
a file called ``fix_esp_data.py`` in the ``code`` directory.

.. code-block:: python

    import os.path


    datadir = os.path.join('..', 'data')
    in_name = os.path.join(datadir, 'jgra20797-sup-0003-ds01.txt')
    out_name = os.path.join(datadir, 'jgra20797-sup-0003-ds01_FIXED.txt')
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
create a file called ``jgra20797-sup-0003-ds01_FIXED.txt`` in the ``data``
directory.

File fixed, we can load and begin examining the data.  Change to the
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

>>> fname = os.path.join(datadir, 'jgra20797-sup-0003-ds01_FIXED.txt')

creates a variable holding the full path to the fixed file.

>>> import numpy

The import statement imports any installed `module <http://docs.python.org/tutorial/modules.html>`_, just as if it were in the standard library. Here we import the very useful :mod:`numpy` module, which is a prerequisite for SpacePy and useful in its own right.

>>> esp_fluxes = numpy.loadtxt(fname, skiprows=14, usecols=[1])

:func:`~numpy.loadtxt` makes it easy to load data from a file into a
numpy :class:`~numpy.ndarray`, a very useful data
container. ``skiprows`` skips the header information, and specifying
only column 1 (first column is column 0) with ``usecols`` will only
load the fluxes for 1.8-3.5MeV. We only load the fluxes at this point
because they can be represented as floats, which numpy arrays store
very efficiently.

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

>>> esp_times = numpy.loadtxt(fname, skiprows=14, usecols=[0],
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
        fname = os.path.join(datadir, 'jgra20797-sup-0003-ds01_FIXED.txt')
        esp_fluxes = numpy.loadtxt(fname, skiprows=14, usecols=[1])
        convert = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d')
        esp_times = numpy.loadtxt(fname, skiprows=14, usecols=[0],
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
OMNI. :func:`~spacepy.omni.get_omni` returns data for a particular
set of times. In this case, we want hourly data, covering 1989 through 
2010 (we'll cut it down to size later). :func:`~spacepy.time.tickrange`
allows us to specify a start time, stop time, and time step.

>>> import spacepy.omni
>>> import spacepy.time
>>> times = spacepy.time.tickrange('1989-01-01', '2011-01-01',
...                                datetime.timedelta(hours=1))
>>> d = spacepy.omni.get_omni(times)
>>> vsw = d['velo']
>>> vsw_times = d['UTC']

We'll also load the esp data:

>>> import common
>>> esp_times, esp_flux = common.load_esp()

Even though we have not installed ``common.py``, the ``import``
statement finds it because it is in the current directory.

``load_esp`` returns a `tuple
<http://docs.python.org/tutorial/datastructures.html#tuples-and-sequences>`_,
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

>>> import scipy.stats
>>> esp_flux_av = numpy.empty(shape=esp_flux.shape, dtype=esp_flux.dtype)
>>> for i in range(len(esp_flux_av)):
...     esp_flux_av[i] = scipy.stats.nanmean(esp_flux[max(i - 13, 0):i + 14])

:func:`numpy.empty` creates an empty array, taking the ``shape`` and
``dtype`` from the ``esp_flux`` array. ``empty`` does not initialize
the data in the array, so it is essentially random junk; use
:func:`~numpy.zeros` to create an array filled with zeros.

:func:`len` returns the length of an array, and :func:`range` then
iterates over each number from 0 to length minus 1, i.e. the entire
array. Each element is then set to a 27-day average: from 13 days
before a day's measurement through 13 days after. (Python slices do
not include the last element listed; they are half-open). Note that
these slices can happily run off the end of the ``esp_flux`` array,
but we use :func:`max` to ensure the first index does not go negative.
(Negative indices have special meaning in Python.)

:func:`~scipy.stats.nanmean` takes the mean of a numpy array,
but skips any elements with a value of "not a number" (nan), which is
often used for fill.  (This is our first exposure to the :mod:`scipy`
module.)

For the solar wind averaging, the times need to cover the 24 * 13.5 = 324
hours previous, and 324 hours following (non-inclusive). There is also a 
more efficient way than using an explicit loop:

>>> vsw_av = numpy.fromiter((scipy.stats.nanmean(vsw[max(0, i - 324):i + 324])
...                         for i in range(len(vsw))),
...			    count=len(vsw), dtype=vsw.dtype)

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

Making Figure 1
===============

To actually plot, we need access to the :mod:`~matplotlib.pyplot` module:

>>> import matplotlib.pyplot as plt
>>> plt.ion()

This alternate form of the import statement shouldn't be overused (it can
make code harder to read by masking the origin of functions), but is
conventional for matplotlib.

:func:`~matplotlib.pyplot.ion` turns on interactive mode so plots appear
and are updated as they're created.

>>> plt.semilogy(esp_times, 10 ** esp_flux_av, 'b')
>>> plt.draw()
>>> plt.draw()

:func:`~matplotlib.pyplot.semilogy` creates a semilog plot, log
on the Y axis. The first two arguments are a list of X and Y values;
after that there are many options to specify formatting (such as the
color, used here.)

The ESP fluxes are stored as the log of the flux; ``**`` is the
exponentiation operator so the (geometric!) average is plotted
properly.

:func:`~matplotlib.pyplot.draw` draws the updated plot; sometimes it
needs to be called repeatedly. Use it whenever you want the plot updated;
it will not be included from here on.

>>> plt.xlabel('Year', weight='bold')
>>> plt.ylabel('Electron Flux\n1.8-3.5 MeV', color='blue', weight='bold')
>>> plt.ylim(1e-2, 10)
(0.01, 10)

:func:`~matplotlib.pyplot.xlabel` and :func:`~matplotlib.pyplot.ylabel`
set the labels for the axes. Note the newline (``\n``) in the string for
the Y label. :func:`~matplotlib.pyplot.ylim` sets the lower and upper
limits for the Y axis; there is, of course, :func:`~matplotlib.pyplot.xlim`
as well.

These are the simplest, although not most flexible, ways to work with plots.
To produce the full Figure 1, we'll move out of interactive mode:

>>> plt.ioff()
>>> plt.show()

:func:`~matplotlib.pyplot.ioff` turns off interactive mode. Once
interactive mode is off, :func:`~matplotlib.pyplot.show` displays
the full plot, including controls for panning, zooming, etc. Until
the plot is closed, nothing further can happen in the Python window.

>>> fig = plt.figure(figsize=[11, 8.5])

:func:`~matplotlib.pyplot.figure` creates a new
:class:`~matplotlib.figure.Figure`; the size specified here is
US-letter paper, landscape orientation.

>>> ax = fig.add_subplot(111)

:meth:`~matplotlib.figure.Figure.add_subplot` creates an
:class:`~matplotlib.axes.Axes` object, which can contain an actual
plot. ``111`` here means that the figure will have 1 subplot and the
new subplot should be in position (1, 1); more on this later.

>>> fluxline = ax.plot(esp_times, 10 ** esp_flux_av, 'b')

:meth:`~matplotlib.axes.Axes.plot` puts the relevant data into the
plot; again specifying a blue line. It returns a list of 
:class:`~matplotlib.lines.Line2D` objects, which we save for later
use.

>>> ax.set_yscale('log')

:meth:`~matplotlib.axes.Axes.set_yscale` switches the Y axis between
log and linear (:meth:`~matplotlib.axes.Axes.set_xscale` for the X axis).

>>> ax.set_ylim(1e-2, 10)
>>> ax.set_xlabel('Year', weight='bold')
>>> ax.set_ylabel('Electron Flux\n1.8-3.5 MeV', color='b', weight='bold')

:meth:`~matplotlib.axes.Axes.set_ylim` (and 
:meth:`~matplotlib.axes.Axes.set_xlim`),
:meth:`~matplotlib.axes.Axes.set_xlabel`, and
:meth:`~matplotlib.axes.Axes.set_ylabel` function much as above, but
operate on a particular :class:`~matplotlib.axes.Axes` object.

>>> ax2 = ax.twinx()

:meth:`~matplotlib.axes.Axes.twinx` establishes a second
Y axis (two values twinned on one X axis) on the same plot.

>>> vswline = ax2.plot(vsw_times, vsw_av, 'r')
>>> ax2.set_ylim(300, 650)
>>> ax2.set_ylabel('Solar Wind Speed', color='r', rotation=270, weight='bold')

The resulting :class:`~matplotlib.axes.Axes` object has all the
methods that we've used before. Note ``rotation`` on
:meth:`~matplotlib.axes.Axes.set_ylabel` to make the text run
top-to-bottom rather than bottom-to-top.

>>> ax.set_xlim(esp_times[0], esp_times[-1])

Since the solar wind data extends beyond the ESP data, this sets
the X axis to match the ESP data. Note ``-1`` to refer to the last
element of the array.

>>> leg = ax.legend([fluxline[0], vswline[0]], ['Flux', 'Vsw'],
...                 loc='upper left', frameon=False)

:meth:`~matplotlib.axes.Axes.legend`, as may be expected, creates a 
:class:`~matplotlib.legend.Legend` on the axes. The first parameter is
a list of the matplotlib objects to make a legend for; since the
plotting commands return these, we can pass them back in. Each plotting
command returns a *list*. In this case we just take the 0th element of
each list since we know there's only one line from each plotting command.
The second parameter is the text used to annotate each line.

>>> fluxtext, vswtext = leg.get_texts()
>>> fluxtext.set_color(fluxline[0].get_color())
>>> vswtext.set_color(vswline[0].get_color())

The default text color is black, so we use
:meth:`~matplotlib.legend.Legend.get_texts` to get the
:class:`~matplotlib.text.Text` objects for the annotations. Again, we
know there are two (we just created the legend). Then
:meth:`~matplotlib.text.Text.set_color` sets the color based on the
the existing color for each line (:meth:`~matplotlib.lines.Line2D.get_color`).

To see the results:

>>> plt.show()

Close the window when done. Now we want to save the output:

>>> fig_fname = os.path.join('..', 'plots', 'fig1a.eps')
>>> fig.savefig(fig_fname)

:meth:`~matplotlib.figure.Figure.savefig` saves the figure, in this case
as an encapsulated PostScript file (to the ``plots`` directory).

Let's tweak a few things. For one, there's a lot of padding around the figure,
which can make it difficult to properly scale for publication. The way around
this is to specify a :class:`~matplotlib.transforms.Bbox` (bounding box),
basically the lower left and upper right corners (in inches) to include
in the saved figure. Getting this right tends to be a matter of trial and error.
(:meth:`~matplotlib.figure.Figure.get_tightbbox` is supposed to help with this,
but it doesn't quite work yet.)

>>> import matplotlib.transforms
>>> bob = matplotlib.transforms.Bbox([[0.52, 0.35], [10.5, 7.95]])
>>> fig.savefig(fig_fname, bbox_inches=bob, pad_inches=0.0)

Better, but all the text is awfully small. Once the figure is fit in the paper
it'll be really small. And the font isn't that great.

>>> import matplotlib
>>> matplotlib.rcParams['axes.unicode_minus'] = False
>>> matplotlib.rcParams['text.usetex']= True
>>> matplotlib.rcParams['font.family'] = 'serif'
>>> matplotlib.rcParams['font.size'] = 14
>>> bob = matplotlib.transforms.Bbox([[0.4, 0.35], [10.7, 7.95]])
>>> fig.savefig(fig_fname, bbox_inches=bob, pad_inches=0.0)

Now the font is bigger and it's rendered using TeX, which should match
the body of the paper better (assuming the paper is in LaTeX). The
larger font means tweaking the bounding box. ``unicode_minus`` fixes a
problem where negative numbers on the axis don't render properly in
TeX. Matplotlib has many more options for `customization
<http://matplotlib.sourceforge.net/users/customizing.html>`_.

The end result is a nice figure that can be printed full-size, put in
a PDF, or included directly in a paper.

Now we need the bottom half of Figure 1. From
`SIDC <http://www.sidc.be/silso/versionarchive>`_, download the "Monthly mean total sunspot number" (``monthssn.dat``). Put it in the ``data``
directory.

>>> import bisect
>>> import datetime
>>> monthfile = os.path.join(common.datadir, 'monthssn.dat')
>>> convert = lambda x: datetime.datetime.strptime(x, '%Y%m')
>>> ssn_data = numpy.genfromtxt(monthfile, skip_header=2400, usecols=[0, 2, 3],
...                             converters={0: convert}, dtype=numpy.object,
...                             skip_footer=24)
>>> idx = bisect.bisect_left(ssn_data[:, 0], datetime.datetime(1989, 1, 1))
>>> ssn_data = ssn_data[idx:]
>>> ssn_times = ssn_data[:, 0]
>>> ssn = numpy.asarray(ssn_data[:, 1], dtype=numpy.float64)
>>> smooth_ssn = numpy.asarray(ssn_data[:, 2], dtype=numpy.float64)
>>> ssn_times += datetime.timedelta(days=15)

Much of this should be familiar. :func:`~numpy.genfromtxt` is a little more
flexible than :func:`~numpy.loadtxt`; here it allows the skipping of lines
at the end as well as the beginning (skipping 200 years at the start, 2 at 
the end, where data are provisional.) Here we load both times and the
sunspot numbers in the same command so that if any lines don't load, they 
will not wind up in any of the arrays.

:mod:`bisect` provides fast functions for searching in sorted data;
:func:`~bisect.bisect_left` is roughly a find-the-position-of function.
Having found the position of the start of 1989, we then keep times
from then on (specifying a start index without a stop index in Python
means "from start to end of the list.") Note that, although ``bisect``
is meant to work on lists, it also works fine on numpy arrays; this is a
common feature of Python known as
`duck typing <http://en.wikipedia.org/wiki/Duck_typing#In_Python>`_.

We then use :func:`~numpy.asarray`
to convert the ``ssn`` and ``smooth_ssn`` columns to float arrays. Note
the slice notation: ``[:, 0]`` means take all indices of the first dimension
(line number) and only the 0th index of the second dimension (column in the
line). Finally, we use :class:`~datetime.timedelta` to shift the date
associated with a month from the beginning to roughly the middle of the month.
Adding a scalar to an array does an element-wise addition.

>>> import matplotlib.figure
>>> fig = plt.figure(figsize=[11, 8.5],
...                  subplotpars=matplotlib.figure.SubplotParams(hspace=0.1))
>>> ax = fig.add_subplot(211)

When creating the figure this time, we use
:class:`~matplotlib.figure.SubplotParams` to choose a slightly smaller
vertical spacing between adjacent subplots. Tweaking ``SubplotParams``
also provides an alternative to tweaking bounding boxes.

Then we create a subplot with the information that there will be 2 rows, 1
column, and this is the first subplot. Now everything acting on ax, above,
can be repeated, although we skip setting the xlabel since only the bottom
axis will be labeled.

>>> fluxline = ax.plot(esp_times, 10 ** esp_flux_av, 'b')
>>> ax.set_yscale('log')
>>> ax.set_ylim(1e-2, 10)
>>> ax.set_ylabel('Electron Flux\n1.8-3.5 MeV', color='b', weight='bold')
>>> ax2 = ax.twinx()
>>> vswline = ax2.plot(vsw_times, vsw_av, 'r')
>>> ax2.set_ylim(300, 650)
>>> ax2.set_ylabel('Solar Wind Speed', color='r', rotation=270, weight='bold')
>>> ax.set_xlim(esp_times[0], esp_times[-1])
>>> leg = ax.legend([fluxline[0], vswline[0]], ['Flux', 'Vsw'],
...                 loc='upper left', frameon=False)
>>> fluxtext, vswtext = leg.get_texts()
>>> fluxtext.set_color(fluxline[0].get_color())
>>> vswtext.set_color(vswline[0].get_color())

Then we move on to adding the solar wind:

>>> ax3 = fig.add_subplot(212, sharex=ax)

This adds another subplot, the second in the 2x1 array. Its x axis is
shared with the existing ``ax``. (This is poorly documented; see this
`example
<http://matplotlib.sourceforge.net/examples/pylab_examples/shared_axis_demo.html>`_)

>>> plt.setp(ax.get_xticklabels(), visible=False)
>>> plt.setp(ax2.get_xticklabels(), visible=False)

:func:`~matplotlib.pyplot.setp` sets a
property. :meth:`~matplotlib.axes.Axes.get_xticklabels` returns all the
tick labels (:class:`~matplotlib.text.Text`) for the x axis; ``setp``
then sets ``visible`` to ``False`` for all of them. This hides the
labeling on the axis for the upper subfigure.

>>> ax3.set_xlabel('Year', weight='bold')
>>> ax3.set_ylabel('Sunspot Number', weight='bold')
>>> smoothline = ax3.plot(ssn_times, smooth_ssn, lw=2.0, color='k')
>>> ssnline = ax3.plot(ssn_times, ssn, color='k', linestyle='dotted')

There is nothing new here except for the specifications of ``linewidth``
and ``linestyle``; see :meth:`~matplotlib.axes.Axes.plot` for details.
Note ``k`` as the abbreviation for black (to avoid confusion with blue.)

>>> leg2 = ax3.legend([ssnline[0], smoothline[0]],
...                   ['Sunspot Number', 'Smoothed SSN'],
...                   loc='upper right', frameon=False)
>>> ax3.set_ylim(0, 200)
>>> ax3.set_xlim(esp_times[0], esp_times[-1])

>>> fig_fname = os.path.join('..', 'plots', 'fig1.eps')
>>> fig.savefig(fig_fname, bbox_inches=bob, pad_inches=0.0)

All of this has been seen for the top half of figure 1.

Following is the complete code to reproduce Figure 1.

.. code-block:: python

    import bisect
    import datetime
    import os.path

    import common
    import matplotlib
    import matplotlib.figure
    import matplotlib.pyplot as plt
    import matplotlib.transforms
    import numpy
    import scipy
    import scipy.stats
    import spacepy.omni
    import spacepy.time


    matplotlib.rcParams['axes.unicode_minus'] = False
    matplotlib.rcParams['text.usetex']= True
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 14
    bob = matplotlib.transforms.Bbox([[0.4, 0.35], [10.7, 7.95]])

    times = spacepy.time.tickrange('1989-01-01', '2011-01-01',
                                   datetime.timedelta(hours=1))
    d = spacepy.omni.get_omni(times)
    vsw = d['velo']
    vsw_times = d['UTC']
    esp_times, esp_flux = common.load_esp()
    esp_flux_av = numpy.empty(shape=esp_flux.shape, dtype=esp_flux.dtype)
    for i in range(len(esp_flux_av)):
        esp_flux_av[i] = scipy.stats.nanmean(esp_flux[max(i - 13, 0):i + 14])
    vsw_av = numpy.fromiter((scipy.stats.nanmean(vsw[max(0, i - 324):i + 324])
                             for i in range(len(vsw))),
                             count=len(vsw), dtype=vsw.dtype)
    monthfile = os.path.join(common.datadir, 'monthssn.dat')
    convert = lambda x: datetime.datetime.strptime(x, '%Y%m')
    ssn_data = numpy.genfromtxt(monthfile, skip_header=2400, usecols=[0, 2, 3],
                                converters={0: convert}, dtype=numpy.object,
                                skip_footer=24)
    idx = bisect.bisect_left(ssn_data[:, 0], datetime.datetime(1989, 1, 1))
    ssn_data = ssn_data[idx:]
    ssn_times = ssn_data[:, 0]
    ssn = numpy.asarray(ssn_data[:, 1], dtype=numpy.float64)
    smooth_ssn = numpy.asarray(ssn_data[:, 2], dtype=numpy.float64)
    ssn_times += datetime.timedelta(days=15)

    fig = plt.figure(figsize=[11, 8.5],
                     subplotpars=matplotlib.figure.SubplotParams(hspace=0.1))
    ax = fig.add_subplot(211)
    fluxline = ax.plot(esp_times, 10 ** esp_flux_av, 'b')
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 10)
    ax.set_ylabel('Electron Flux\n1.8-3.5 MeV', color='b', weight='bold')
    ax2 = ax.twinx()
    vswline = ax2.plot(vsw_times, vsw_av, 'r')
    ax2.set_ylim(300, 650)
    ax2.set_ylabel('Solar Wind Speed', color='r', rotation=270, weight='bold')
    ax.set_xlim(esp_times[0], esp_times[-1])
    leg = ax.legend([fluxline[0], vswline[0]], ['Flux', 'Vsw'],
                    loc='upper left', frameon=False)
    fluxtext, vswtext = leg.get_texts()
    fluxtext.set_color(fluxline[0].get_color())
    vswtext.set_color(vswline[0].get_color())

    ax3 = fig.add_subplot(212, sharex=ax)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax3.set_xlabel('Year', weight='bold')
    ax3.set_ylabel('Sunspot Number', weight='bold')
    smoothline = ax3.plot(ssn_times, smooth_ssn, lw=2.0, color='k')
    ssnline = ax3.plot(ssn_times, ssn, color='k', linestyle='dotted')
    leg2 = ax3.legend([ssnline[0], smoothline[0]],
                      ['Sunspot Number', 'Smoothed SSN'],
                      loc='upper right', frameon=False)
    ax3.set_ylim(0, 200)
    ax3.set_xlim(esp_times[0], esp_times[-1])

    fig_fname = os.path.join('..', 'plots', 'fig1.eps')
    fig.savefig(fig_fname, bbox_inches=bob, pad_inches=0.0)

.. _appendix:

Appendix: Fixing the ESP data file
==================================
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

>>> in_name = os.path.join(datadir, 'jgra20797-sup-0003-ds01.txt')
>>> out_name = os.path.join(datadir, 'jgra20797-sup-0003-ds01_FIXED.txt')
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

`pop <http://docs.python.org/tutorial/datastructures.html#more-on-lists>`_ returns one element from a list, and deletes it from the list. Using ``0`` pops off the first element, and :meth:`~file.write` writes a string to a file. ``+`` can be used to concatenate two strings together. Since :meth:`~str.split` removed the newlines, they need to be readded.

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
<http://docs.python.org/tutorial/controlflow.html#if-statements>`_
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
