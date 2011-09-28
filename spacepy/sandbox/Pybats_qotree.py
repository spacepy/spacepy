#!/usr/bin/env python

'''
Let's try to make a OQ tree.  QO tree.  Whatever.

Original by Dan Welling, cleanup and subclassing by Brian Larsen

'''

import numpy as np
from qotree import *


def pybats_decision(tree=None, branchNum=None):
    a = np.sqrt(tree[branchNum].npts)
    if (int(a) == a) and (a > 0 ):
        dx = (np.max(grid[0,:][tree[branchNum].locs]) \
                  - np.min(grid[0,:][tree[branchNum].locs])) \
                  / (np.sqrt(tree[branchNum].npts)-1)
        if not np.isnan(dx):
            if int(np.log2(dx))==np.log2(dx):
                # An NxN block can be considered a "leaf" or stopping point
                # if Dx is an integer power of 2.  Leafs must "know" the
                # indices of the x,y points located inside of them as a
                # grid and also know the bounding coords of the grid cells.
                #tree[branchNum].isLeaf = True
                tree[branchNum].locs=self[i].locs.reshape( (a, a) )
                tree[branchNum].dx = dx
                tree[branchNum].cells = np.meshgrid(
                    np.arange(tree[branchNum].lim[0], tree[branchNum].lim[1]+dx, dx),
                    np.arange(tree[branchNum].lim[2], tree[branchNum].lim[3]+dx, dx))
                return False
        else:
            return False
    else:
        return True
