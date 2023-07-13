######################
#                    #
# Michele Martinazzo #
#                    #
######################

# 

try: import numpy as np
except ImportError as msg:  raise SystemExit (str(msg) + '\nimport numpy (numeric python) failed!')

def _optics_geometry(a,xi,phi0,theta0):
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

    phi    = phi0   + (dphi)*180/np.pi
    theta  = theta0 - (dtheta)*180/np.pi

    return phi, theta

####################################################################################################################################

def optics_shape(phi0, theta0, shape = "circular", a = 0, xi = 0):

    '''
    '''
    
    if shape == "circular":
        a  = a/2.
        xi = np.linspace(0, 2*np.pi,1001)

    elif shape == "custom":
        a  = np.asarray(a)
        xi = np.asarray(xi)

    else:
        raise SystemExit ('\nError: the allowed shape are circular or custom!') 

    phi, theta = _optics_geometry(a,xi,phi0,theta0)

    return phi, theta

####################################################################################################################################

if __name__ == '__main__':
    pass



