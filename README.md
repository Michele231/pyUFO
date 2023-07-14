![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/logo.png) 
# pyUFO

pyUFO is a Python library specifically designed to assist with the management of satellite field of view (FOV) data. The current version of the library provides the ability to project a general FOV shape onto a sphere (using a specific function), returning the latitude and longitude of the FOV contour. Future versions will aim to expand its capabilities, allowing for the projection of FOV shapes onto ellipsoids and geoids. 

## Installation

To install pyUFO library it is sufficient to clone this repository and run the setup script. 

```bash
git clone https://github.com/Michele231/pyUFO.git
cd pyUFO
pip install .
```

## Usage

pyUFO allows to use the function "on_sphere()" to create a projection on the sphere surface.

on_sphere() INPUTS:
1. ssp_lat  : latitude of the sub satellite point 
2. ssp_lon  : longitude of the sub satellite point 
3. hsat     : Height of the satellite (km)
4. phi0     : FOV central optical axis azimuth angle (with respect to N) (°)
5. theta0   : FOV central optical axis zenith angle (90° = nadir) (°)
6. r_opt    : Opening radius of the optics (mrad) (np array)
7. xi_opt   : Angle of the optics (°) (np array)
8. shape    : Shape of the optics (circular or custom)
9. r_sphere : Radius of the sphere (km)

on_sphere() OUTPUTS:
1. latf  : latitude of the FOV contour on the sphere
2. lonf  : longitude of the FOV contour on the sphere

### Example:

```python
from pyUFO import on_sphere

# fancy FOV shape (in polar coordinates)
xi_opt = np.linspace(1, 360,360)
r_opt = 150*(np.sin((xi*np.pi/180)*4))**2+150

# compute the lat lon coordinates of the FOV
lat, lon = on_sphere(ssp_lat = -45, ssp_lon = 0, hsat = 800, 
                    phi0 = 0, theta0 = 45, r_opt=r_opt, xi_opt = xi_opt, shape='custom')
```

### Some Figures:

FOV shape             |  on_sphere()
:-------------------------:|:-------------------------:
Custom shape | theta0 = 45°
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/fov_geometry_custom.png)  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/45d_obs_custom.png)
Circular shape | theta0 = 90° (Nadir)
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/fov_geometry.png)  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/nadir_obs1.png)
Circular shape | theta0 = 30°
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/fov_geometry.png "Circular shape")  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/30d_obs1.png "theta0 = 30°")


## Contributing

Pull requests are welcome. If you find any bug or error, please let me know!

## License

[MIT](https://choosealicense.com/licenses/mit/)
