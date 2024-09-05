The radbelt Module
==================

The radiation belt module currently includes a simple radial
diffusion code as a class. Import the module and instatiate a radbelt object::

>>> import spacepy.radbelt as sprb
>>> rb = sprb.RBmodel()

Add a time grid for a particular period that you are interested in::

>>> rb.setup_ticks('2002-02-01T00:00:00', '2002-02-10T00:00:00', 0.25)

This will automatically lookup required geomagnetic/solar wind conditions
for that period. Run the diffusion solver for that setup and plot the
results::

>>> rb.evolve()
>>> rb.plot()


The Data Assimilation Module
============================

This module includes data assimilation capabilities, through the
assimilation class. The class assimilates data for the radiation belt model
using the Ensemble Kalman Filter. The algorithm used is the SVD method
presented by Evensen in 2003 (Evensen, G., Ocean dynamics, 53, pp.343--367,
2003). To compensate for model errors, three inflation algorithms are
implemented. The inflation methodology is specified by the inflation
argument, where the options are the following:

   * inflation = 0: Add model error (perturbation for the ensemble) around
     model state values only where observations are available (DEFAULT).

   * inflation = 1: Add model error (perturbation for the ensemble) around
     observation values only where observations are available.

   * inflation = 2: Inflate around ensemble average for EnKF.

Prior to assimilation, a set of data values has to be specified by setting the
start and end dates, and time step, using the ``setup_ticks`` function of the
radiation belt model::

>>> import spacepy
>>> import datetime
>>> from spacepy import radbelt

>>> start = datetime.datetime(2002,10,23)
>>> end = datetime.datetime(2002,11,4)
>>> delta = datetime.timedelta(hours=0.5)
>>> rmod.setup_ticks(start, end, delta, dtype='UTC')

Once the dates and time step are specified, the data is added using the
``add_PSD`` function (NOTE: This requires a database available from the SpacePy team)::

>>> rmod.add_PSD()

The observations are averaged over the time windows, whose interval is give by
the time step. Once the dates and data are set, the assimilation is performed
using the ``assimilate`` function::

>>> rmod.assimilate(inflation=1)

This function will add the PSDa values, which are the analysis state of
the radiation belt using the observations within the dates. To plot the
analysis simply use the ``plot`` function::

>>> rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

Additionally, to create a summary plot of the observations use the ``plot_obs``
function within the radbelt module. For reference, the last closed drift shell,
Dst, and Kp are all included. These can be disabled individually using the
corresponding Boolean kwargs.

The clims kwarg can be used to manually set the color bar range.  To use, set
it equal to a two-element list containing minimum and maximum log :sub:`10` value to
plot.  Default action is to use [0,10] as the log :sub:`10` of the color range.  This
is good enough for most applications. The title of the top most plot defaults
to 'Summary Plot' but can be customized using the title kwarg.

The figure object and all three axis objects (PSD axis, Dst axis, and Kp axis)
are all returned to allow the user to further customize the plots as necessary.
If any of the plots are excluded, None is returned in their stead.

Example::

>>> rmod.plot_obs(clims=[-10,-6],Lmax=False,Kp=False,Dst=False,title='Observations Plot')

This command would create the summary plot with a color bar range of 10 :sup:`-10`
to 10 :sup:`-6`.  The Lmax line, Kp and Dst values would be excluded.  The title of
the topmost plot (phase space density) would be set to 'Observations Plot'.
