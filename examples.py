import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# call the library
from pyUFO import on_sphere, on_ellipsoid

# circular
#lat, lon = on_sphere(45, 0, 800, 180, 90, r_opt=240)

# fancy FOV shape
xi = np.linspace(1, 360,3600)
r_opt = 15 #mrad

# compute the lat lon coordinates of the FOV
lats, lons = on_sphere(ssp_lat = 0.0, ssp_lon = 8, hsat = 800, 
                    phi0 = 0, theta0 = 60, r_opt=r_opt, xi_opt = xi, shape='custom')

# compute the lat lon coordinates of the FOV
late, lone = on_ellipsoid(ssp_lat = 0.0, ssp_lon = 8, hsat = 800, 
                    phi0 = 0, theta0 = 60, r_opt=r_opt, xi_opt = xi, shape='custom')


# plot on map
plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(1, figsize=(7,7),
                      subplot_kw=dict(projection=ccrs.PlateCarree()))

ax.set_extent([7.2, 9.2, 3.2, 5])
ax.coastlines(resolution='10m')
ax.set_title("")

# Create a scatter plot
scatplot = ax.scatter(lons, lats, c='blue', s=1.5,
                      transform=ccrs.PlateCarree())
scatplot = ax.scatter(lone, late, c='red', s=1.5,
                      transform=ccrs.PlateCarree())



# Sort out gridlines and their density
#xticks_extent = list(np.arange(-50, 50, 10))
#yticks_extent = list(np.arange(-20, 80, 5))
xticks_extent = list(np.arange(7.2, 9.2, 0.2))
yticks_extent = list(np.arange(0, 10, 0.2))

gl = ax.gridlines(linewidths=0.1)
gl.xlabels_top = False
gl.xlabels_bottom = True
gl.ylabels_left = True
gl.ylabels_right = False
gl.xlocator = mticker.FixedLocator(xticks_extent)
gl.ylocator = mticker.FixedLocator(yticks_extent)
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

plt.show()