######################
#                    #
# Michele Martinazzo #
#                    #
######################

# This set of functions solve the cone-sphere intersection problem.

import numpy as np

def shape_to_viewing_angle(a,xi,phi0,theta0):
    '''
    This function transform the shape of the optic in a set line of sight.

    INPUTS:
            - a     :
            - xi    :
            - phi0  : central optical axis azimuth angle (with respect to N) (°)
            - theta0: central optical axis zenith angle (90 = nadir) (°)
    OUTPUTS:
            - phi   :
            - theta :

    '''
    ca     = np.cos(a)
    sa     = np.sin(a)

    b      = sa*np.cos(xi)
    c      = sa*np.sin(xi)

    dphi   = np.arctan(b/ca)
    dtheta = np.arctan(c)

    phi    = phi0 + (dphi)*180/np.pi
    theta  = theta0 - (dtheta)*180/np.pi

    return phi,theta



