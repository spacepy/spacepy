#!/usr/bin/env python
'''
Tools for creating interactive PyBats GUIs.
'''

class ClickTracer(object):
    '''
    ClickTracers, given an axis containing a 2D plot and the pybats.bats.Bats2d
    object used to create it, allows the user to click to add a streamline 
    passing through the point clicked.
    '''

    def __init__(self, ax, bats, xfield, yfield, debug=False):
        # Check arguments.
        assert type(ax).__module__ == 'matplotlib.axes', \
            'First argument must be a matplotlib.axes.Axes (sub)object.'
        assert type(bats).__name__ == 'Bats2d', \
            '2nd argument must be a Bats2d object.'
        self.bats   = bats
        self.ax     = ax
        self.xfield = xfield
        self.yfield = yfield
        self.debug  = debug
        self.connected = False
        self.connect()

    def _do_trace(self, event, **kwargs):
        if event.inaxes != self.ax: 
            return
        if self.debug:
            print('click at %.3f, %.3f' % (event.xdata, event.ydata))
        stream = self.bats.get_stream(event.xdata, event.ydata,
                                      self.xfield, self.yfield)
        ylims = self.ax.get_ylim()
        xlims = self.ax.get_xlim()
        self.nowline = stream.plot(self.ax, alpha=0.25)
        self.ax.set_ylim(ylims)
        self.ax.set_xlim(xlims)
        self.ax.figure.canvas.draw()

    def connect(self):
        '''
        Activate trace-on-click.
        '''
        if self.connected:
            return
        self.cid = self.ax.figure.canvas.mpl_connect(
            'button_press_event',self._do_trace)
        #self.kid = self
        self.connected = True

    def disconnect(self):
        '''
        Deactivate trace-on-click.
        '''
        if not self.connected:
            return
        self.ax.figure.canvas.mpl_disconnect(self.cid)
        self.connected = False
