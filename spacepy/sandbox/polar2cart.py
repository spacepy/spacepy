

"""
Have a look at this and see if we want to inclide it somewhere

Brian Larsen
30Aug2011
"""

import numpy as np

def polar2cart(r, theta):
    """
    Convert polar coordinates, r, theta to Cartesian x, y

    This is done in a few ways to be as fast as it can be

    Parameters
    ==========
    r : array-like float
        radius values
    theta: array-like float
        angle values in radians

    Returns
    =======
    out : array-like
        the x, y coordinates of the input

    Examples
    ========
    r = [3.]*1000000
    theta = [np.pi/4.]*1000000
    start = time.time()
    x, y = polar2cart(r, theta)

    See Also
    ========
    cart2polar
    """
    try:
        import numexpr as ne
    except ImportError:
        x = np.multiply(r, np.cos(theta))
        y = np.multiply(r, np.sin(theta))
        return x, y
    x = ne.evaluate("r * cos(theta)")
    y = ne.evaluate("r * sin(theta)")
    return x, y


def cart2polar(x, y):
    """
    Convert Cartesian x, y to polar coordinates, r, theta

    This is done in a few ways to be as fast as it can be

    Parameters
    ==========
    x : array-like float
        x values
    y: array-like float
        y values

    Returns
    =======
    out : array-like
        the r, theta coordinates of the input in radians

    Examples
    ========
    x = np.zeros(1000000)+3.
    y = np.zeros(1000000)+4.
    start = time.time()
    r, theta = cart2polar(x, y)

    See Also
    ========
    polar2cart
    """
    try:
        import numexpr as ne
    except ImportError:
        r = np.sqrt(np.power(x, 2) + np.power(y, 2))
        theta = np.arctan2(y, x)
        return r, theta
    r = ne.evaluate("sqrt(x**2 + y**2)")
    theta = ne.evaluate("arctan2(y,x)")
    return r, theta


def polar2cart2(r, theta):
    x = np.multiply(r, np.cos(theta))
    y = np.multiply(r, np.sin(theta))
    return x, y

def cart2polar2(x, y):
    r = np.sqrt(np.power(x, 2) + np.power(y, 2))
    theta = np.arctan2(y, x)
    return r, theta




if __name__ == "__main__":
    import numexpr as ne
    import time
    r = [3.]*1000000
    theta = [np.pi/4.]*1000000
    start = time.time()
    x, y = polar2cart(r, theta)
    end = time.time()
    print("%4.2f: 1e6 lists polar2cart numexpr" % (end-start) )

    r = [3.]*1000000
    theta = [np.pi/4.]*1000000
    start = time.time()
    x, y = polar2cart2(r, theta)
    end = time.time()
    print("%4.2f: 1e6 lists polar2cart numpy" % (end-start) )

    ##########################
    r = np.zeros(1000000)+3
    theta = np.zeros(1000000)+ np.pi/4.
    start = time.time()
    x, y = polar2cart(r, theta)
    end = time.time()
    print("%4.2f: 1e6 arrays polar2cart numexpr" % (end-start) )

    r = np.zeros(1000000)+3
    theta = np.zeros(1000000)+ np.pi/4.
    start = time.time()
    x, y = polar2cart2(r, theta)
    end = time.time()
    print("%4.2f: 1e6 arrays polar2cart numpy" % (end-start) )

    ##########################

    x = [3]*1000000
    y = [4]*1000000
    start = time.time()
    r, theta = cart2polar(x, y)
    end = time.time()
    print("%4.2f: 1e6 lists cart2polar numexpr" % (end-start) )

    x = [3]*1000000
    y = [4]*1000000
    start = time.time()
    r, theta = cart2polar2(x, y)
    end = time.time()
    print("%4.2f: 1e6 lists cart2polar numpy" % (end-start) )

    ##########################

    x = np.zeros(1000000)+3.
    y = np.zeros(1000000)+4.
    start = time.time()
    r, theta = cart2polar(x, y)
    end = time.time()
    print("%4.2f: 1e6 arrays cart2polar numexpr" % (end-start) )

    x = np.zeros(1000000)+3.
    y = np.zeros(1000000)+4.
    start = time.time()
    r, theta = cart2polar2(x, y)
    end = time.time()
    print("%4.2f: 1e6 arrays cart2polar numpy" % (end-start) )
