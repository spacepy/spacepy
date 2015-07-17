#!/usr/bin/env python

'''
A class for building QO trees for :class:`~spacepy.pybats.bats.Bats2d` 
objects.
'''

class QTree(object):
    '''
    Base class for Quad/Oct tree objects assuming cell-centered grid points
    in a rectangular non-regular layout.
    '''
    
    def __init__(self, grid):
        '''
        Build QO Tree for input grid.  Grid should be a NxM numpy array where
        N is the number of dimensions and M is the number of points.
        '''
        from numpy import sqrt, where, arange, lexsort
        (self.d, self.npoints) = grid.shape

        if(self.d != 2):
            raise NotImplementedError("Sorry, QTrees are for 2D grids only.")
        self.tree = {}
        
        # Find limits of cell-centered grid.
        # Get values and locations of grid max/min
        minloc, xmin = grid[0,:].argmin(), grid[0,:].min()
        maxloc, xmax = grid[0,:].argmax(), grid[0,:].max()
        # Distance from min X value is grid size at that point.
        r = sqrt((grid[0,:]-xmin)**2. + (grid[1,:]-grid[1,:][minloc])**2.)
        ml=(where(r>0, r, 10000)).min()
        # Actual min is not cell center, but cell boundary.
        xmin-=ml/2.
        # Repeat for Xmax.
        r = sqrt((grid[0,:]-xmax)**2. + 
                 (grid[1,:]-grid[1,:][maxloc])**2.)
        ml=(where(r>0, r, 10000)).min()
        xmax+=ml/2.
        
        # Repeat for Ymin/max.
        minloc, ymin = grid[1,:].argmin(),grid[1,:].min()
        maxloc, ymax = grid[1,:].argmax(),grid[1,:].max()

        r = sqrt((grid[1,:]-ymin)**2. + (grid[0,:]-grid[0,:][minloc])**2.)
        ml=(where(r>0, r, 10000)).min(); ymin-=ml/2.
        r = sqrt((grid[1,:]-ymax)**2. + (grid[0,:]-grid[0,:][maxloc])**2.)
        ml=(where(r>0, r, 10000)).min(); ymax+=ml/2.

        # Some things all trees should know about themselves.
        self.nleafs=0
        self.aspect_ratio = (xmax-xmin)/(ymax-ymin)

        # Use spatial range of grid to seed root of QTree.
        self[1]=[xmin,xmax,ymin,ymax]
        self.locs = lexsort( (grid[0,:], grid[1,:]) )
        self._spawn_kids(grid)
        self.nbranch=len(list(self.keys()))

    def _spawn_kids(self, grid, i=1):
        '''
        Internal recursive method for populating tree.
        '''
        from numpy import sqrt, max, min, log2, arange, meshgrid
        self[i].locs=self.locs[(grid[0,:][self.locs]>self[i].lim[0]) &
                               (grid[0,:][self.locs]<self[i].lim[1]) &
                               (grid[1,:][self.locs]>self[i].lim[2]) &
                               (grid[1,:][self.locs]<self[i].lim[3]) ]
        self[i].npts = self[i].locs.size
        # Bats blocks are almost always 8x8; 
        # blocks must have no less than 64 points.
        # Stop branching as soon as possible (e.g. combine blocks of like dx).
        a = log2(self[i].npts/64.0) # where a is npts=2^a * 64
        if (int(a) == a): 
            # integer 'a' implies correct number of points to be a "leaf".
            # Approximate dx assuming a proper block.
            xmax=max(grid[0,:][self[i].locs])
            xmin=min(grid[0,:][self[i].locs])
            dx = (xmax-xmin) / (sqrt(self[i].npts)-1)
            # Count points along x=xmax and x= xmin.  These are equal in Leafs. 
            nxmin=len(grid[0,:][self[i].locs][grid[0,:][self[i].locs]==xmin])
            nxmax=len(grid[0,:][self[i].locs][grid[0,:][self[i].locs]==xmax])
            # Define leaf as area of constant dx (using approx above) 
            # or npts=64 (min level.)
            if (a==0) or (nxmax==nxmin==sqrt(self[i].npts)):
                # An NxN block can be considered a "leaf" or stopping point
                # if above criteria are met.  Leafs must "know" the 
                # indices of the x,y points located inside of them as a 
                # grid and also know the bounding coords of the grid cells.
                self[i].isLeaf = True
                a=sqrt(self[i].npts)
                self[i].locs=self[i].locs.reshape( (a, a) )
                self[i].dx = dx
                self[i].cells = meshgrid(
                    arange(self[i].lim[0], self[i].lim[1]+dx, dx),
                    arange(self[i].lim[2], self[i].lim[3]+dx, dx))
                self.nleafs+=1
                return

        # If above criteria are not met, this block is 
        # not a constant-resolution zone.
        # Subdivide section into four new ones (8 if oct tree)
        dx = (self[i].lim[1] - self[i].lim[0])/2.0
        x=[self[i].lim[0], self[i].lim[0]+dx, self[i].lim[0]+dx, self[i].lim[0]]
        y=[self[i].lim[2], self[i].lim[2], self[i].lim[2]+dx, self[i].lim[2]+dx]
        for j, k in enumerate(range(self.ld(i), self.rd(i)+1)):
            self[k] = [x[j],x[j]+dx,y[j],y[j]+dx]
            self._spawn_kids(grid,k)


    def find_leaf(self, x, y, i=1):
        '''
        Recursively search for and return the index of the leaf that 
        contains the input point x, y.
        '''
        l=self[i].lim
        # if point is in this block...
        if(l[0]<=x<=l[1])and(l[2]<=y<=l[3]):
            # ...and it's a leaf, return this block.
            if self[i].isLeaf: return i
            # ...and it's a branch, dig deeper.
            else:
                for j in range(self.ld(i), self.rd(i)+1):
                    answer = self.find_leaf(x,y,i=j)
                    if answer: return answer
        else:
            return False

    def __getitem__(self, key):
        return self.tree[key]

    def __setitem__(self, key, value):
        self.tree[key] = Branch(value)

    def __contains__(self, key):
        return key in self.tree

    def __iter__(self):
        return iter(self.tree)

    def keys(self):
        return list(self.tree.keys())

    def mom(self, k):
        return (k+2**self.d-2)/2**self.d
    def leftdaughter(self, k):
        return k*2**self.d-2**self.d+2
    def rightdaughter(self,k):
        return k*2**self.d+1
    # convenience:
    rd = rightdaughter
    ld = leftdaughter

    def plot_res(self, ax, DoLabel=True, tagLeafs=False):
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
        for key in self.tree:
            if self[key].isLeaf:
                if self[key].dx in res_colors:
                    dx_vals[self[key].dx] = 1.0
                    color=res_colors[self[key].dx]
                else:
                    color='k'
                self[key].plot_res(ax, fc=color, label=key*tagLeafs)

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
    def __init__(self, lim):
        '''
        lim should be a 4 element list of the
        dimensional boundaries of the branch.
        '''
        self.isLeaf = False
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

    def plot_res(self, ax, fc='gray', label=False):
        if not self.isLeaf: return
        from matplotlib.patches import Polygon
        from numpy import array
        l=self.lim
        verts = array([
                [l[0],l[2] ], [ l[1],l[2]],
                [l[1],l[3] ], [ l[0],l[3]]])
        poly = Polygon(verts, True, ec=None, fc=fc, lw=0.0001)
        if label:
            x = l[0]+(l[1]-l[0])/2.0
            y = l[2]+(l[3]-l[2])/2.0
            ax.text(x, y, label) 
        ax.add_patch(poly)
