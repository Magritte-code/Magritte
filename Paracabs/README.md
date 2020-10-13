# Paracabs

[![Build Status](https://travis-ci.com/FredDeCeuster/Paracabs.svg?branch=master)](https://travis-ci.com/FredDeCeuster/Paracabs)

Parallelization and acceleration abstractions for code and performance portability.
By abstracting away the implementation specifics, we provide a limited but unified interface for parallelization and acceleration.


## Abstraction back ends

### Multi-threading
- OpenMP

### Message passing
- MPI

### Acceleration
- Cuda
- Sycl


_This library was built for use in the 3D radiative transfer code Magritte._
_It is largly inspired by, but not as complete as:_
- _Grid, by Peter Boyle et al._
- _Hemi, by Mark Harris et al._
- _Eigen_
