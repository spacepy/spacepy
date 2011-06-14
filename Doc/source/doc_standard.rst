Documentation Standard
======================

SpacePy uses Sphinx (http://sphinx.pocoo.org/) as its documentation system

On top of Sphinx SpacePy uses the following extensions:
 * 'sphinx.ext.autodoc'
 * 'sphinx.ext.doctest''
 * 'sphinx.ext.intersphinx'
 * 'sphinx.ext.todo'
 * 'sphinx.ext.pngmath'
 * 'sphinx.ext.ifconfig'
 * 'sphinx.ext.viewcode'
 * 'numpydoc'
 * 'sphinx.ext.inheritance_diagram'
 * 'sphinx.ext.autosummary'
 * 'sphinx.ext.extlinks'


So what do I need to do in my code?
-----------------------------------
Since we are using the 'numpydoc' extension there are fixed heading that may
appear in your documentation block, there are a few things to note:
 * No other headings can appear there
 * I believe that no reStructuredText commands can either (e.g. .. Note:)

Allowed headings
~~~~~~~~~~~~~~~~
**Always use**
 * Parameters
 * Returns

**Use as needed**
 * Attributes
 * Raises
 * Warns
 * Other Parameters
 * See Also
 * Notes
 * Warnings
 * References
 * Examples
 * Methods

**No need to use**
 * Summary
 * Extended Summary
 * index

**Do not use**
 * Signature

**Examples**
 * Use them, but they must be fully stand alone, the user can just copy-paste
   exactly what is there and it should work (think, and probably run doctest on them)

Method Example
--------------
This code from toolbox shows what a method should look like in your code

    .. code-block:: python

        def logspace(min, max, num, **kwargs):
            """
            Returns log spaced bins.  Same as numpy logspace except the min and max are the ,min and max
            not log10(min) and log10(max)

            Parameters
            ==========
            min : float
                minimum value
            max : float
                maximum value
            num : integer
                number of log spaced bins

            Other Parameters
            ================
            kwargs : dict
                additional keywords passed into matplotlib.dates.num2date

            Returns
            =======
            out : array
                log spaced bins from min to max in a numpy array

            Notes
            =====
            This function works on both numbers and datetime objects

            Examples
            ========
            >>> import spacepy.toolbox as tb
            >>> tb.logspace(1, 100, 5)
            array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])
            """
            from numpy import logspace, log10
            if isinstance(min, datetime.datetime):
                from matplotlib.dates import date2num, num2date
                return num2date(logspace(log10(date2num(min)), log10(date2num(max)), num, **kwargs))
            else:
                return logspace(log10(min), log10(max), num, **kwargs)


Which then renders as:

    .. autofunction:: toolbox.logspace
