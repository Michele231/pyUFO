######################
#                    #
# Michele Martinazzo #
#                    #
######################

# 

try: import numpy as np
except ImportError as msg:  raise SystemExit (str(msg) + '\nimport numpy (numeric python) failed!')

from .optic import optics_shape

def _sat2earth(x0,y0,z0,phi,theta):
    '''
    This function tranform the viewing direction from the
    satellite basis to the Earth basis.
    '''
    if (z0 == 0): z0 = z0 + 1e-9

    # satellite basis
    v0 = np.array([x0, y0, z0])                                # up
    v0 = v0/np.linalg.norm(v0)                                 # Frobenius Norm
    if np.abs(z0) <= 1e-7:
        v1 = np.sign(x0)*np.array([0, 0, -1])                             
    else:
        v1 = np.array([1, y0/x0, 
        -(x0/z0 + (y0*y0)/(x0*z0))])                           # N
    v1 = np.sign(x0*z0)*v1/np.linalg.norm(v1)
    v2 = np.cross(v0,v1)                                       # E

    # view direction in the satellite basis
    d_view = np.array([-np.sin(theta), np.cos(theta)*np.cos(phi+np.pi), 
                    np.cos(theta)*np.sin(phi+np.pi)])

    # transformation matrix: satellite --> earth
    T            = np.ndarray((3,3),buffer=np.asarray([v0,v1,v2])).T
    d_earth_view = T.dot(d_view)

    return d_earth_view

####################################################################################################################################

def _find_intersection_with_x(d_earth_view, r_sphere, x0,y0,z0):
    '''
    This
    '''
    
    # Solving the line-sphere intersection
    a = np.squeeze(d_earth_view[0,:])
    a[np.abs(a)<=1e-8] = np.sign(a[np.abs(a)<=1e-8])*1e-8
    b = np.squeeze(d_earth_view[1,:])
    b[np.abs(b)<=1e-8] = np.sign(b[np.abs(b)<=1e-8])*1e-8
    c = np.squeeze(d_earth_view[2,:])
    c[np.abs(c)<=1e-8] = np.sign(c[np.abs(c)<=1e-8])*1e-8

    k1 = (b/a)*x0-y0
    k2 = (c/a)*x0-z0

    alph = 1 + (b/a)**2 + (c/a)**2
    beta = -2.*(k1*b/a+k2*c/a)
    gamm = k1**2+k2**2-r_sphere**2

    # find the real solutions
    Delt = (beta/alph)**2-4*(gamm/alph)
    alph[Delt<0] = np.nan
    beta[Delt<0] = np.nan
    gamm[Delt<0] = np.nan
    Delt[Delt<0] = np.nan

    x1 = -0.5*((beta/alph)+np.sqrt(Delt))
    x2 = -0.5*((beta/alph)-np.sqrt(Delt))
    y1 = (b/a)*x1 - (b/a)*x0 + y0
    y2 = (b/a)*x2 - (b/a)*x0 + y0
    z1 = (c/a)*x1 - (c/a)*x0 + z0
    z2 = (c/a)*x2 - (c/a)*x0 + z0

    return x1,x2,y1,y2,z1,z2

####################################################################################################################################

def _x2latlon(xf,yf,zf,r_sphere):
    '''
    This
    '''
    
    thetaf = np.arccos(zf/r_sphere)
    phif   = np.sign(yf)*np.arccos(xf/np.sqrt(xf**2+yf**2))

    latf   = 90 - 180*thetaf/np.pi
    lonf   = 180*phif/np.pi

    return latf, lonf

####################################################################################################################################

def fov_on_sphere(ssp_lat, ssp_lon, hsat, 
                    phi0, theta0, r_sphere = 6371):
    '''
    This
    '''

    phi0   = np.pi*phi0/180 
    theta0 = np.pi*theta0/180

    # checks
    if (ssp_lat <= -90 or ssp_lat >= 90):
        raise SystemExit ('\nError: lat of the ssp has to be >-90 and <90')
    if (ssp_lon < -180 or ssp_lat > 180):
        raise SystemExit ('\nError: lon of the ssp has to be >=-180 and <=180')

    #np.asarray(phi0)
    #np.asarray(theta0)

    # computation of x0, y0 and z0 (position of the satellite)
    lat_ang = np.pi/2 - np.pi*ssp_lat/180
    lon_ang = np.pi*ssp_lon/180

    R_tot   = hsat + r_sphere
    x0      = R_tot*np.sin(lat_ang)*np.cos(lon_ang)
    y0      = R_tot*np.sin(lat_ang)*np.sin(lon_ang)
    z0      = R_tot*np.cos(lat_ang)

    # satellite view vector
    d_earth_view = _sat2earth(x0,y0,z0,phi0,theta0)

    # Solving the line-sphere intersection
    x1,x2,y1,y2,z1,z2 = _find_intersection_with_x(d_earth_view, 
                                            r_sphere,x0,y0,z0)

    # find the closest set of points to x0,y0,z0
    d1 = np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
    d2 = np.sqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)

    xf,yf,zf = x2,y2,z2
    xf[d2>d1],yf[d2>d1],zf[d2>d1] = x1[d2>d1],y1[d2>d1],z1[d2>d1]


    # transformation to lat-lon coordinates 
    latf, lonf = _x2latlon(xf,yf,zf,r_sphere)

    return latf, lonf

####################################################################################################################################
    
def on_sphere(ssp_lat, ssp_lon, hsat,
              azimuth, zenith, r_opt, xi_opt = 0, shape = "circular",
              r_sphere = 6371):

    '''
    FOV projection on a sphere

    INPUTS:
            - ssp_lat  : latitude of the sub satellite point 
            - ssp_lon  : longitude of the sub satellite point 
            - hsat     : Height of the satellite (km)
            - azimuth  : central optical axis azimuth angle (with respect to N) (°)
            - zenith   : central optical axis zenith angle (90 = nadir) (°)
            - r_opt    : Opening radius of the optics (mrad) (np array)
            - xi_opt   : Angle of the optics (°) (np array)
            - shape    : Shape of the optics (circular or custom)
            - r_sphere : Radius of the sphere (km)
    OUTPUTS:
            - latf     : np array containing the latitudes of the fov
            - lonf     : np array containing the longitudes of the fov

    '''


    # checks
    if (zenith < 0 or zenith > 90):
        raise SystemExit ('\nError: zenith has to be >=0 and <=90')
    if (azimuth < -180 or azimuth > 180):
        raise SystemExit ('\nError: azimuth has to be >=-180 and <=180')

    phi, theta = optics_shape(azimuth,zenith,shape,r_opt,xi_opt)
    latf, lonf = fov_on_sphere(ssp_lat, ssp_lon, hsat, 
                    phi, theta, r_sphere)

    return latf, lonf

####################################################################################################################################

if __name__ == "__main__":
    pass
