#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot: SpacePy specialized plotting routines

This package provides classes to plot various types of space physics data.

Authors: Brian Larsen and Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov


Copyright 2011 Los Alamos National Security, LLC.

"""
from matplotlib.patches import Wedge
import matplotlib.pyplot as plt

__all__ = ["spectrogram", "utils"]



def dual_half_circle(center=(0,0), radius=1.0,
                     sun_direction='right', ax=None, colors=('w','k'),
                     **kwargs):
    """
    Plot two half circles to a plot with the specified face colors and
    rotation. This is normal to use to denote the sun direction in
    magnetospheric science plots.

    Other Parameters
    ----------
    center : array-like, 2 elements
        Center in data coordinates of the circles, default (0,0)
    radius : float
        Radius of the circles, defualt 1.0
    sun_direction : string or float
        The rotation direction of the first (white) circle. Options are ['down',
        'down right', 'right', 'up left', 'up right', 'up', 'down left', 'left']
        or an angle in degrees counter-clockwise from up. Default right.
    ax : matplotlib.axes
        Axis to plot the circles on.
    colors : array-like, 2 elements
        The two colors for the circle fill. The First number is the light and
        second is the dark.
    **kwargs : other keywords
        Other keywords to pass to matplotlib.patches.Wedge

    Returns
    -------
    out : tuple
        Tuple of the two wedge objects

    Examples
    --------
    >>> import spacepy.plot
    >>> spacepy.plot.dual_half_circle()

    """
    sun_dict = {'left':90,
                'right':-90,
                'up':0,
                'down':180,
                'up right':-45,
                'up left':45,
                'down right':45+180,
                'down left':-45+180,
                }
    try:
        angle = sun_dict[sun_direction]
    except KeyError:
        try:
            angle = float(sun_direction)
        except ValueError:
            raise(ValueError("sun_direction was not understood, must be a float or {0}".format(sun_dict.keys())))
                
    theta1, theta2 = angle, angle + 180
    if ax is None:
        ax = plt.gca()
        
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return (w1, w2)

