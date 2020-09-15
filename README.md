# Paracabs
Parallelization and acceleration abstractions for code and performance portability.
By abstracting away the implementation specifics, we provide a limited but unified interface for parallelization and acceleration.

_This library was initially built for use in the 3D radiative transfer code Magritte._

Based on:
- Grid
- Nvidia's/Mark Harris' hemi
- Eigen

## Available abstraction back ends

Multi-threading
- OpenMP

Message passing
- MPI

Acceleration
- Cuda
- Sycl
