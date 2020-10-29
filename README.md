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
- CUDA
- SYCL

## Disclaimer
This library was built for use in the 3D radiative transfer code Magritte.

## References
This library is largely inspired by, but not as complete as:
- Grid, by Peter Boyle et al.
- Hemi, by Mark Harris et al.
- Eigen
