import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# call the library
from pyUFO import on_sphere, on_ellipsoid

# circular
#lat, lon = on_sphere(45, 0, 800, 180, 90, r_opt=240)

# fancy circular FOV with Gaussian dispersion
xi = np.linspace(0, 360,361)
rr    = np.linspace(0, 200,401)
sigma = 14

# plot on map
plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(1, figsize=(7,7),
                      subplot_kw=dict(projection=ccrs.PlateCarree()))

ax.set_extent([6.6, 9.8, 38.8, 41.6])
ax.coastlines(resolution='10m')
ax.set_title("")

# Create a scatter plot
#scatplot = ax.scatter(lons, lats, c='blue', s=1.5,
#                      transform=ccrs.PlateCarree())
for r in rr:
	color      = (0.2, 0.5, 0.5,0.5*np.exp(-(r**2/(2*sigma**2))))
	late, lone = on_ellipsoid(ssp_lat = 35.0, ssp_lon = 15, hsat = 800, 
                    azimuth = 45, zenith = 45, r_opt=r, xi_opt = xi, shape='custom',
                    semimajor_ax = 6378.137, semiminor_ax = 6356.752)
	scatplot   = ax.scatter(lone, late, c=color, s=1.5,
	                      transform=ccrs.PlateCarree())
late, lone = on_ellipsoid(ssp_lat = 35.0, ssp_lon = 15, hsat = 800, 
                azimuth = 45, zenith = 45, r_opt=14, xi_opt = xi, shape='custom',
                semimajor_ax = 6378.137, semiminor_ax = 6356.752)
scatplot   = ax.scatter(lone, late, c=(0.1, 0.4, 0.4), s=1.5,
	                      transform=ccrs.PlateCarree())


# Sort out gridlines and their density
#xticks_extent = list(np.arange(-50, 50, 10))
#yticks_extent = list(np.arange(-20, 80, 5))
xticks_extent = list(np.arange(6.6, 9.8, 0.2))
yticks_extent = list(np.arange(38, 42, 0.2))

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

fig, ax = plt.subplots(1, figsize=(7,7))
ax.plot(xi,late,c='blue')
#ax.plot(xi,lone,c='red')
plt.show()

