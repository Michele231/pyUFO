import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# call the library
from pyUFO import on_sphere

# circular
#lat, lon = on_sphere(45, 0, 800, 180, 90, r_opt=240)

# fancy FOV shape
xi = np.linspace(1, 360,360)
r_opt = 150*(np.sin((xi*np.pi/180)*4))**2+150

# compute the lat lon coordinates of the FOV
lat, lon = on_sphere(ssp_lat = -45, ssp_lon = 0, hsat = 800, 
                    phi0 = 0, theta0 = 45, r_opt=r_opt, xi_opt = xi, shape='custom')


# plot on map
plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(1, figsize=(7,7),
                      subplot_kw=dict(projection=ccrs.PlateCarree()))

ax.set_extent([-50, 50, -70, 0.0])
ax.coastlines(resolution='10m')
ax.set_title("")

# Create a scatter plot
scatplot = ax.scatter(lon, lat, c='blue', s=1.5,
                      transform=ccrs.PlateCarree())



# Sort out gridlines and their density
xticks_extent = list(np.arange(-50, 50, 10))
yticks_extent = list(np.arange(-70, 0, 5))

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