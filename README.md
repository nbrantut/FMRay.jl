# FMRay.jl

This module is mostly a port of some key elements of the C++ programmes FaATSO, which implements:
 - Fast Marching Method to compute arrival times in 3D, VTI (elliptical) wave velocity structures,
 - Ray tracing and computation of derivatives of arrival times,
 - Event location method by gridsearch methods with least L1 norm.

Immediate applications can be written to solve the Local Earthquake Tomography problem (as was done originally in FaATSO), but this is left for future work, and maybe another package.

Since this package is mostly a port to Julia, the key original paper describing this code is:

Brantut, N. (2018), Time-resolved tomography using acoustic emissions in the laboratory, and application to snadstone compaction, *Geophys. J. Int.*, 213, 2177-2192.

and references therein.

---

The package is not registered. You should add it using:
` pkg> add https://github.com/nbrantut/FMRay.jl`

Then test using:
` pkg> test`

---

## Velocity model

FMRay main tool is an implementation of the Fast Marching Method to solve the Eikonal equation in 3D, possibly anisotropic media.

The type of anisotropy is restricted to VTI, with the following parameterisation:
$$V(\theta) = V_v \cos^2(\theta) + V_h \sin^2(\theta),$$
where $V_v$ and $V_h$ are the vertical and horizontal phase velocities, respectively, $\theta$ is the phase angle, and $V(\theta)$ is the phase velocity. Note that the solver is not valid for anisotropy that is too large.

## Usage

A velocity model is defined in a `Grid` structure, which contains the grid spacing (only regular grids can be used), vertical and horizontal phase velocities, and the coordinates of the origin. The velocity structure are 3-dimensional arrays.

The origin point must be entered as the coordinates of that point in the default coordinate system, where the point indexed (1,1,1) is at position (0.0, 0.0, 0.0).

An example with 2D uniform velocity equal to 1.0, and grid spacing 0.1:

````julia
using FMRay

G = Grid(0.1, fill(1.0, (100,100,1)))

T = march(CartesianIndex(1,1,1), G)
````

The function `march` is the Eikonal solver. The first argument is the index of the source, and must be entered as a CartesianIndex.
