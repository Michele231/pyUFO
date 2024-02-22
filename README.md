![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/logo.png) 

<p align="center">
  <img title="v0.1" alt="v0.1" src="https://img.shields.io/badge/version-v0.0.1-informational?style=flat-square">
  <img title="MIT License" alt="license" src="https://img.shields.io/badge/license-MIT-informational?style=flat-square">
	<img title="python" alt="python3.11" src="https://img.shields.io/badge/python-3.11-informational?style=flat-square"><br/>
	<img title="Code size" alt="code size" src="https://img.shields.io/github/languages/code-size/Michele231/pyUFO?color=red">
	<img title="Repo size" alt="repo size" src="https://img.shields.io/github/repo-size/Michele231/pyUFO?color=red">
</p>

***

# pyUFO

pyUFO is a Python library specifically designed to assist with the management of satellite field of view (FOV). The current version of the library provides the ability to project a general FOV shape onto a sphere or onto a ellipsoid (analytically), returning the latitude and longitude of the FOV contour. 

## Installation

To install pyUFO library it is sufficient to clone this repository and run the setup script. 

```bash
git clone https://github.com/Michele231/pyUFO.git
cd pyUFO
pip install .
```

## Usage

Two functions are currently available in pyUFO: "on_sphere()" and "on_ellipsoid()".

### on_sphere()

**INPUTS**:
1. ssp_lat  : Latitude of the sub-satellite point. 
2. ssp_lon  : Longitude of the sub-satellite point.
3. hsat     : Height of the satellite perpendicular to the surface (km).
4. phi0     : Observation direction (central optical axis) azimuth angle (counterclockwise with respect to N) (°).
5. theta0   : Observation direction (central optical axis) zenith angle (90° = nadir) (°).
6. r_opt    : Polar coordinate radius for the description of the optic shape (mrad) (np array).
7. xi_opt   : Polar coordinate angle for the description of the optic shape (°) (np array).
8. shape    : Shape of the optics (circular or custom).
9. r_sphere : Radius of the sphere (km).

**OUTPUTS**:
1. latf  : latitude of the FOV contour on the sphere.
2. lonf  : longitude of the FOV contour on the sphere.

### on_ellipsoid()

**INPUTS**:
1. ssp_lat      : Latitude of the sub-satellite point. 
2. ssp_lon      : Longitude of the sub-satellite point.
3. hsat         : Height of the satellite perpendicular to the surface (km).
4. phi0         : Observation direction (central optical axis) azimuth angle (counterclockwise with respect to N) (°).
5. theta0       : Observation direction (central optical axis) zenith angle (90° = nadir) (°).
6. r_opt        : Polar coordinate radius for the description of the optic shape (mrad) (np array).
7. xi_opt       : Polar coordinate angle for the description of the optic shape (°) (np array).
8. shape        : Shape of the optics (circular or custom).
9. semimajor_ax : Equatorial radius (semi-major axis) (km).
10. semiminor_ax : Polar radius (semi-minor axis) (km).

**OUTPUTS**:
1. latf  : latitude of the FOV contour on the ellispoid.
2. lonf  : longitude of the FOV contour on the ellispoid.

Note that, for ellipsoids, the conversion of ECEF coordinates to latitude involves a circular relationship, and it is solved iteratively using Ferrari's solution. This method may encounter numerical instability when the difference between the semi-major and semi-minor axes becomes significant. To circumvent this issue, it is recommended to maintain the relative difference between these axes below 15%. However the ECEF coordinates are always well determinated, and can be retrieved from the main function (xf,yf,zf).

### Example:

```python
from pyUFO import on_sphere

# fancy instrument optics shape (in polar coordinates)
xi = np.linspace(1, 360,361)
r  = 150*(np.sin((xi*np.pi/180)*4))**2+150 # mrad

# compute the lat lon coordinates of the FOV on the Earth sphere approximation
lat, lon = on_sphere(ssp_lat = -45, ssp_lon = 0, hsat = 800, 
                    phi0 = 0, theta0 = 45, r_opt = r, xi_opt = xi, shape='custom')

# compute the lat lon coordinates of the FOV on the Earth IERS (2003) reference ellipsoid
lat, lon = on_ellipsoid(ssp_lat = -45, ssp_lon = 0, hsat = 800, 
                    phi0 = 0, theta0 = 45, r_opt = r, xi_opt = xi, shape='custom')
```

### Some Figures:

FOV shape             |  on_sphere()
:-------------------------:|:-------------------------:
Custom shape | theta0 = 45°
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/fov_geometry_custom.png)  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/45d_obs_custom.png)
Circular shape | theta0 = 30°
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/fov_geometry.png "Circular shape")  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/30d_obs1.png "theta0 = 30°")

Circular FOV          |  on_ellipsoid()
:-------------------------:|:-------------------------:
Gaussian Dispersion | Latitude vs Polar coordinate angle xi_opt
![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/gaussian_dispersion.png)  |  ![alt text](https://github.com/Michele231/pyUFO/blob/main/figures/lat_vs_xi.png)




## Contributing

Pull requests are welcome. If you find any bug or error, please let me know!

## License

[MIT](https://choosealicense.com/licenses/mit/)

***

Authors: Michele Martinazzo.
