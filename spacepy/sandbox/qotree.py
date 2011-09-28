#!/usr/bin/env python

'''
Let's try to make a OQ tree.  QO tree.  Whatever.

Original by Dan Welling, cleanup and subclassing by Brian Larsen

'''

import numpy as np

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



def leftdaughter(dim, k):
    return k*2**dim-2**dim+2

def rightdaughter(dim, k):
    return k*2**dim+1

def mother(dim, k):
    return (k+2**dim-2)/2**dim

def boxes_through_level(dim, level):
    return ((2**dim)**level - 1)/(2**dim - 1)


class QTree(dict):
    '''
    Base class for Quad/Oct tree objects assuming cell-centered grid points
    in a rectangular non-regular layout.
    '''

    def __init__(self, grid, grid_edges=False, max_depth = 3, decision_func = None):
        '''
        Build QO Tree for input grid.  Grid should be a NxM numpy array where
        N is the number of dimensions and M is the number of points.

        Parameters
        ==========
        grid, numpy.ndarray
            Grid should be a NxM numpy array where N is the number of dimensions and M is the number of points.
        '''
        if (grid.shape[0] != 2):
            raise(NotImplementedError("Sorry, QTrees are for 2D grids only."))

        self.decision_func = decision_func

        self.max_depth = max_depth
        self.grid = grid

        (self.d, self.npoints) = grid.shape

        if not grid_edges:
            #####
            # toolbox.bin_center_to_edges can be used to make centers into edges
            #####
            # Find limits of cell-centered grid.
            # Get values and locations of grid max/min
            minloc, xmin = grid[0,:].argmin(), grid[0,:].min()
            maxloc, xmax = grid[0,:].argmax(), grid[0,:].max()
            # Distance from min X value is grid size at that point.
            r = np.sqrt((grid[0,:]-xmin)**2. + (grid[1,:]-grid[1,:][minloc])**2.)
            ml=r[r>0].min()
            # Actual min is not cell center, but cell boundary.
            xmin -= ml/2.
            # Repeat for Xmax.
            r = np.sqrt((grid[0,:]-xmax)**2. + (grid[1,:]-grid[1,:][maxloc])**2.)
            ml=r[r>0].min()
            xmax += ml/2.

            # Repeat for Ymin/max.
            minloc, ymin = grid[1,:].argmin(), grid[1,:].min()
            maxloc, ymax = grid[1,:].argmax(), grid[1,:].max()

            r = np.sqrt((grid[1,:]-ymin)**2. + (grid[0,:]-grid[0,:][minloc])**2.)
            ml=r[r>0].min()
            ymin -= ml/2.
            r = np.sqrt((grid[1,:]-ymax)**2. + (grid[0,:]-grid[0,:][maxloc])**2.)
            ml=r[r>0].min()
            ymax += ml/2.

            self.aspect_ratio = (xmax-xmin)/(ymax-ymin)

            # Use spatial range of grid to seed root of QTree.
            self[1] = Branch(np.asarray([xmin,xmax,ymin,ymax]))
            self.locs = np.lexsort( (grid[1,:], grid[0,:]) ) # lexsort works from the right
            self._spawn_daughters()
        else:
            raise(NotImplementedError("Sorry, Only allowed to use grid centers so far"))

    def _spawn_daughters(self, i=1):
        '''
        Internal recursive method for populating tree.
        '''
        self[i].locs = self.locs[(self.grid[0,:][self.locs]>self[i].lim[0]) &
                               (self.grid[0,:][self.locs]<self[i].lim[1]) &
                               (self.grid[1,:][self.locs]>self[i].lim[2]) &
                               (self.grid[1,:][self.locs]<self[i].lim[3]) ]
        self[i].npts = self[i].locs.size
        if self.decision_func != None:
            if not self.decision_func(self, i):
                self[i].isLeaf = True
                return
        # Subdivide section into four new ones (8 if oct tree)
        dx = (self[i].lim[1] - self[i].lim[0])/2.0
        x=[self[i].lim[0], self[i].lim[0]+dx, self[i].lim[0]+dx, self[i].lim[0]]
        y=[self[i].lim[2], self[i].lim[2], self[i].lim[2]+dx, self[i].lim[2]+dx]
        for j, k in enumerate(range(self.getleftdaughter(i), self.getrightdaughter(i)+1)):
            self[k] = Branch(np.asarray([x[j],x[j]+dx,y[j],y[j]+dx]))
            # if we are at max depth we don't want to split again
            # this is tested by seeing if k is in the lowest level
            if k <= self.getboxes_through_level(self.max_depth-1):
                self._spawn_daughters(k)
            else:
                self[k].isLeaf = True


    def getboxes_through_level(self, level):
        return boxes_through_level(self.d, level)

    def getmother(self, k):
        return mother(self.d, k)

    def getleftdaughter(self, k):
        return leftdaughter(self.d, k)

    def getrightdaughter(self, k):
        return rightdaughter(self.d, k)

    def plot_res(self, ax, DoLabel=True):
        res_colors={
            1./32.: 'black',
            1./16.: 'darkred',
            1./8. : 'red',
            1./4. : 'orange',
            1./2. : 'yellow',
            1.    : 'green',
            2.    : 'darkblue',
            4.    : 'blue',
            8.    : 'lightblue',
            16.   : 'grey',
            32.   : 'black'}

        dx_vals = {}
        for key in self:
            if self[key].isLeaf:
                self[key].plot_res(ax, fc=res_colors[self[key].dx])
                dx_vals[self[key].dx] = 1.0

        if DoLabel:
            ax.annotate('Resolution:', [1.01,0.99], xycoords='axes fraction',
                        color='k',size='medium')
            for i,key in enumerate(sorted(dx_vals.keys())):
                if key<1:
                    label = '1/%i' % (key**-1)
                else:
                    label = '%i' % key
                ax.annotate('%s $R_{E}$'%label, [1.01,0.95-i*0.05],
                            xycoords='axes fraction', color=res_colors[key],
                            size='x-large')

class Branch(object):
    '''
    Base class for branches/leafs along a QO tree.
    '''
    def __init__(self, lim, isLeaf=False):
        '''
        lim should be a 4 element list of the
        dimensional boundaries of the branch.
        '''
        try:
            assert(len(lim) == 4)
        except AssertionError:
            raise(ValueError("Limits can only be a 4 element array"))
        self.isLeaf = isLeaf
        self.lim = lim

    def plotbox(self, ax, lc='k', **kwargs):
        '''
        Plot a box encompassing the branch lim onto
        axis 'ax'.
        '''
        from matplotlib.collections import LineCollection
        from numpy import array
        l=self.lim
        segs = (
            array([ [l[0],l[2] ], [ l[1],l[2]] ]),
            array([ [l[0],l[3] ], [ l[1],l[3]] ]),
            array([ [l[0],l[2] ], [ l[0],l[3]] ]),
            array([ [l[1],l[2] ], [ l[1],l[3]] ]))
        coll=LineCollection(segs, colors=lc, **kwargs)
        ax.add_collection(coll)
        #ax.plot(l[0:2], [l[2],l[2]], **kwargs)
        #ax.plot(l[0:2], [l[3],l[3]], **kwargs)
        #ax.plot([l[0],l[0]], l[2:],  **kwargs)
        #ax.plot([l[1],l[1]], l[2:],  **kwargs)

    def plot_res(self, ax, fc='gray'):
        if not self.isLeaf: return
        from matplotlib.patches import Polygon
        from numpy import array
        l=self.lim
        verts = array([
                [l[0],l[2] ], [ l[1],l[2]],
                [l[1],l[3] ], [ l[0],l[3]]])
        poly = Polygon(verts, True, ec=None, fc=fc, lw=0.0001)
        ax.add_patch(poly)
