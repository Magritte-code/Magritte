<img src="docs/src/images/Magritte_logo.svg" alt="logo" width="350"/>

[![Build Status](https://travis-ci.com/FredDeCeuster/Magritte.svg?token=j3NNTbFLxGaJNsSoKgCz&branch=master)](https://travis-ci.com/FredDeCeuster/Magritte)
[![Documentation Status](https://readthedocs.com/projects/magritte-magritte/badge/?version=latest&token=7d7b32178f0e9cfbbf7e0967233b680ccf113dd73b8d4b567f3587445c936036)](https://magritte-magritte.readthedocs-hosted.com/en/latest/?badge=latest)
---

_*Preparing for first release...*_

This is the repository for Magritte: a modern open- source software library for
radiation transport simulation. Currently, the emphasis is on 3D radiative
transfer modelling for astrophysics and cosmology, but the techniques could also
be applied more general. Magritte uses a deterministic ray-tracer with a
formal solver that currently focusses on line radiative transfer. See this paper
([arXiv](https://arxiv.org/pdf/1912.08445.pdf),
[MNRAS](https://doi.org/10.1093/mnras/stz3557)) for more details.

Magritte can either be used as a C++ library or as a Python package.


## Installation
First, download the dependencies and configure Magritte using
```bash
bash dependencies/get_dependencies.sh
```
Then, create an anaconda environment from the environment file with
```bash
conda env create -f dependencies/conda_env.yml
```
Afterwards, activate the environment you just created with
```bash
conda activate magritte
```
Now Magritte can be build using
```bash
bash build.sh
```
This builds the library in the `/bin` folder. Note that if you try to build
Magritte from outside the `magritte_env` conda environment compilation might
fail or the generated library might not work in python due to version mismatches.
Therefore as a general rule: **always build and use Magritte from within the
magritte_env conda environment**. To use Magritte's python interface, you should
include the package folder in your python path e.g. by including
```python
from sys import path
path.append("path/to/Magritte")
```
Please have a look at the `build.sh` script for further configuration options.


## Issues
Please report any issues [here](https://github.com/UCL/Magritte/issues).


## Dependencies
* `CMake`, for building;
* `Eigen`, for linear algebra;
* `Paracabs`, for parallelisation and acceleration;
* `pybind11`, for interfacing with python;
* `googletest`, for testing.


## Papers
* _Magritte II: Adaptive ray-tracing, mesh construction and reduction_
([arXiv](https://arxiv.org/abs/2011.14998), [MNRAS](https://doi.org/10.1093/mnras/staa3199));
* _Magritte I: Non-LTE atomic and molecular line modelling_
([arXiv](https://arxiv.org/abs/1912.08445),
[MNRAS](https://doi.org/10.1093/mnras/stz3557)).


## Cite
Please contact the authors of the papers referenced above if you want to use
Magritte in your research. We are currently working on documentation and
examples to facilitate its independent use. Until then, please contact us.
