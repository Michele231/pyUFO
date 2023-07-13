######################
#                    #
# Michele Martinazzo #
#                    #
######################

# 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

try: import numpy as np
except ImportError as msg:  raise SystemExit (str(msg) + '\nimport numpy (numeric python) failed!')

from optic import optics_shape

def _sat2earth(x0,y0,z0,phi,theta):
    '''
    This
    '''


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
    
    phif   = np.arccos(zf/r_sphere)
    thetaf = np.sign(yf)*np.arccos(xf/np.sqrt(xf**2+yf**2))

    latf   = 90 - 180*phif/np.pi
    lonf   = 180*thetaf/np.pi

    return latf, lonf

####################################################################################################################################

def fov_on_sphere(ssp_lat, ssp_lon, hsat, 
                    phi0, theta0, r_sphere = 6371):
    '''
    This
    '''

    # checks
    if (ssp_lat <= -90 or ssp_lat >= 90):
        raise SystemExit ('\nError: lat of the ssp has to be >-90 and <90')
    if (ssp_lon < -180 or ssp_lat > 180):
        raise SystemExit ('\nError: lon of the ssp has to be >=-180 and <=180')

    np.asarray(phi0)
    np.asarray(theta0)

    # computetation of x0, y0 and z0 (position of the satellite)
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



    # just a figure
    #fig = plt.figure()
    # sphere
    #u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
    #x = r_sphere*np.cos(u)*np.sin(v)
    #y = r_sphere*np.sin(u)*np.sin(v)
    #z = r_sphere*np.cos(v)
    #ax = fig.add_subplot(1,1,1, projection='3d')
    #plot = ax.plot_surface(
	#	x, y, z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
#		linewidth=0, antialiased=False, alpha=0.5)
    # sat_point
    #plot = ax.scatter(x0,y0,z0, marker='*')
    # obs direction 1
    #plot = ax.quiver(x0,y0,z0, 5*hsat*d_earth_view[0,0], 
    #               5*hsat*d_earth_view[1,0], 5*hsat*d_earth_view[2,0])
    # obs direction 2
    #plot = ax.quiver(x0,y0,z0, 5*hsat*d_earth_view[0,1], 
    #               5*hsat*d_earth_view[1,1], 5*hsat*d_earth_view[2,1])
    #kkkk=0
    # p1
    #plot = ax.scatter(xf[kkkk],yf[kkkk],zf[kkkk], marker='.')


    #plt.show()

    return latf, lonf

####################################################################################################################################
    
if __name__ == "__main__":

    a,b=optics_shape(1, 1, shape = "circular")
    fov_on_sphere(-45,34,800,np.array([0,np.pi]),np.array([np.pi/3,np.pi/3]))
    #fov_on_sphere(-89,-180,800,np.array([0,]),np.array([0,]))
    #fov_on_sphere(-89,-180,800,np.array([0,]),np.array([0,]))
    #fov_on_sphere(-89,-180,800,[0,2],[0,2])
