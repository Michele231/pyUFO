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
    '''
    a      = a/1000

    ca     = np.cos(a)
    sa     = np.sin(a)

    b      = sa*np.cos(xi)
    c      = sa*np.sin(xi)

    dphi   = np.arctan(b/ca)
    dtheta = np.arcsin(c)

    phi    = phi0   + (dphi)*180/np.pi
    theta  = theta0 - (dtheta)*180/np.pi


    return phi, theta

####################################################################################################################################

def _x2spherical(xf,yf,zf,r_sphere):
    '''
    (90 = nadir)
    '''
    
    thetaf = np.arccos(zf/r_sphere)
    phif   = np.sign(yf)*np.arccos(xf/np.sqrt(xf**2+yf**2))

    theta   = 180*thetaf/np.pi-90
    phi     = 180*phif/np.pi

    return theta, phi

####################################################################################################################################

def optics_shape(phi0, theta0, shape = "circular", a = 0, xi = 0):

    '''
        This function transform the shape of the optic in a set line of sight.

    INPUTS:
            - phi0  : central optical axis azimuth angle (with respect to N) (°)
            - theta0: central optical axis zenith angle (90 = nadir) (°)
            - a     :
            - xi    :
    OUTPUTS:
            - phi   :
            - theta :
    '''
    
    if shape == "circular":
        xi = np.linspace(0, 2*np.pi,1001)

    elif shape == "custom":
        a  = np.asarray(a)
        xi = np.asarray(xi)
        xi = np.pi*xi/180

    else:
        raise SystemExit ('\nError: the allowed shape are circular or custom!') 

    # transform to rad
    phi0   = np.pi*phi0/180
    theta0 = np.pi*theta0/180 # it would be + np.pi/2 but then is removed

    # computing the shape of the optics
    phi, theta = _optics_geometry(a,xi,0,90)
    phi        = np.pi*phi/180
    theta      = np.pi*theta/180
    
    # computing the cartesian reference of frame
    x          = np.sin(theta)*np.cos(phi)
    y          = np.sin(theta)*np.sin(phi)
    z          = np.cos(theta)

    # rotation matrices 
    ry0 = np.array([np.cos(theta0),0,-np.sin(theta0)])
    ry1 = np.array([0,1,0])
    ry2 = np.array([np.sin(theta0),0,np.cos(theta0)])
    Ry  = np.ndarray((3,3),buffer=np.asarray([ry0,ry1,ry2])).T
    rz0 = np.array([np.cos(phi0),np.sin(phi0),0])
    rz1 = np.array([-np.sin(phi0),np.cos(phi0),0])
    rz2 = np.array([0,0,1])
    Rz  = np.ndarray((3,3),buffer=np.asarray([rz0,rz1,rz2])).T

    # rotation
    xyz     = np.ndarray((3,len(x)),buffer=np.asarray([x,y,z]))
    xyz_rot = Rz.dot(Ry.dot(xyz))

    # back to spherical coordinates
    theta_rot, phi_rot = _x2spherical(xyz_rot[0,:],xyz_rot[1,:],xyz_rot[2,:],1)

    return phi_rot, theta_rot

####################################################################################################################################

if __name__ == '__main__':
    pass



