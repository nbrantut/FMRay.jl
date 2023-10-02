# FaVTI.jl

This module is mostly a port of some key elements of the C++ programmes FaATSO, which implements:
 - Fast Marching Method to compute arrival times in 3D, VTI (elliptical) wave velocity structures,
 - Ray tracing and computation of derivatives of arrival times,
 - Event location method by gridsearch methods with least L1 norm.

Immediate applications can be written to solve the Local Earthquake Tomography problem (as was done originally in FaATSO), but this is left for future work, and maybe another package.

Since this package is mostly a port to Julia, the key original paper describing this code is:

Brantut, N. (2018), Time-resolved tomography using acoustic emissions in the laboratory, and application to snadstone compaction, /Geophys. J. Int./, 213, 2177-2192.

and references therein.

---

The package is not registered. You should add it using:
` pkg> add https://github.com/nbrantut/FaVTI.jl`
