
import numpy as np
def PlaneFit(XYZ):
    """
    modified from http://omar.toomuchcookies.net/node/2010/10/residual-of-planefitting-using-svd/

    retuns the orthogonal distance regression plane that minimizes the
    perpendicular distances to the plane for a set of points

    The equation of the a plane is a*x + b*y + c*z + d = 0

    Parameters
    ==========
    XYZ : np.array
        data points in shape nx3

    Returns
    =======
    B : np.array
        coefficients of the plane fit : a*x + b*y + c*z + d = 0

    Examples
    ========
    import numpy as np
    points = np.asarray([[1,-6, 0], [-4, 2, -5], [-2, 4, 1]])
    a, b, c, d = PlaneFit(points) # or just ans = to get an array out
    np.allclose(a * (-2) + b * 4  + c * 1 + d, 0)
    # True

    """

    XYZ = XYZ.transpose()
    rows,npts = XYZ.shape # this is an 3xn array, need to transpose what we have
    if not rows == 3:
            print XYZ.shape
            raise ('data is not 3D')
            return None
    if npts <3:
            raise ('too few points to fit plane')
            return None
    # Set up constraint equations of the form  AB = 0,
    # where B is a column vector of the plane coefficients
    # in the form   b(1)*X + b(2)*Y +b(3)*Z + b(4) = 0.
    t = XYZ.T
    p = (np.ones((npts,1)))
    A = np.hstack([t,p])
    if npts == 3:        # Pad A with zeros
            A = np.vstack((A, np.zeros((1,4))))
    u, d, v = np.linalg.svd(A)  # Singular value decomposition.
    #print v[3,:]
    B = v[3,:]       # Solution is last column of v.
    nn = np.linalg.norm(B[0:3])
    B = B / nn
    return B

