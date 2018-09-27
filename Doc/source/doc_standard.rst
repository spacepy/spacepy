Documentation Standard
======================

SpacePy aims to be a high quality product, and as such we (the SpacePy Team) encourage 
a a high degree of uniformity in the documentation across included modules.  If you are 
contributing to SpacePy, or hope to, please take the time to make your code compliant 
with the documentation standard.

SpacePy uses Sphinx_ to generate its documentation. This allows most of the documentation 
to be built from docstrings in the code, with additional information being provided in
reStructured Text files. This allows easy generation of high-quality, searchable HTML 
documentation.

In addition to Sphinx, SpacePy uses the following extensions:
 * 'sphinx.ext.autodoc'
 * 'sphinx.ext.doctest''
 * 'sphinx.ext.intersphinx'
 * 'sphinx.ext.todo'
 * 'sphinx.ext.imgmath' (falls back to 'sphinx.ext.pngmath' if imgmath is not available)
 * 'sphinx.ext.ifconfig'
 * 'sphinx.ext.viewcode'
 * 'numpydoc'
 * 'sphinx.ext.inheritance_diagram'
 * 'sphinx.ext.autosummary'
 * 'sphinx.ext.extlinks'

.. _Sphinx: http://sphinx.pocoo.org/

So what do I need to do in my code?
-----------------------------------
Since we are using the 'numpydoc' extension there are fixed headings that may
appear in your documentation block. There are a few things to note:
* No other headings can appear in your docstrings
* Most reStructuredText commands cannot appear in your docstrings either (e.g. .. Note:)
* Since 'numpydoc' is not well documented, the best way of finding out what you can do in your docstrings is to look at the source for the SpacePy documentation or the numpy documentation.

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
 * Use them, but they must be fully stand alone; the user should be able to type the exact 
   code in the example and it should work as shown (doctest can help with this)

Function Example
----------------
This code from toolbox shows what a function should look like in your code

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

    .. autofunction:: spacepy.toolbox.logspace
        :noindex:


