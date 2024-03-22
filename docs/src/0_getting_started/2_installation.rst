.. _link-installation:

Installation
############

.. note::

    This is the comprehensive installation guide. For a quick intro, see our
    :ref:`quickstart <link-quickstart>` guide.


Download
********

Magritte has to be compiled from its source code, which can be cloned using:

.. code-block:: shell

    git clone --recursive https://github.com/Magritte-code/Magritte.git

from our `GitHub <https://github.com/Magritte-code/Magritte>`_ repository. This
creates the directory, :literal:`Magritte`, which will be refered to as the
Magritte root directory. Ensure to include the :literal:`--recursive` to also
clone the required submodules.


Dependencies
************

Magritte has several dependencies, some of which are optional.


**Submodules**

* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, version :literal:`3.3.7` or later, for some of the linear algebra;

* `pybind11 <https://github.com/pybind/pybind11>`_, version :literal:`2.2.4` or later, for binding C++ to python;

* `googletest <https://github.com/google/googletest>`_ version :literal:`1.10.0` or later, for some of the tests;

* `Paracabs <https://github.com/Magritte-code/Paracabs>`_, our custom parallelization and acceleration abstractions.

**Required**

* `GCC <https://gcc.gnu.org/>`_, version :literal:`5.0.0` or later, to compile the C++ part of Magritte;
* `CMake <https://cmake.org/>`_, version :literal:`3.18.0` or later, for building the library, organising compilation and linking;


**Optional**

* `Anaconda <https://www.anaconda.com/blog/individual-edition-2020-11>`_, for managing the required Python packages;
* `Gmsh <https://gmsh.info/>`_, version :literal:`4.6.0` or later, for meshing model geometries.


Please note that :literal:`Paracabs` might have further dependencies depending
on which paralellization and acceleration libraries are used. See
:ref:`advanced compilation <link-advanced_compilation>` for further details.


Python packages & Environment
*****************************

To make Magritte useful, some python packages have to be installed which
is best done using a `conda <https://www.anaconda.com/products/individual>`_ environment.

The following python packages are used in Magritte, mainly for io and to create
the model files.

* :mod:`numpy`, to bind the Magritte data;
* :mod:`h5py`, to read and write HDF5 data files;
* :mod:`scipy`, for interpolation and spatial functions such as nearest neighbour calculations;
* :mod:`healpy`, to sample directions from a discretized unit sphere;
* :mod:`astropy`, for unit conversions and physical constants;
* :mod:`meshio`, for reading and writing several types of mesh data structures;
* :mod:`vtk`, for reading and writing vtk files;
* :mod:`pyyaml`, for reading and writing yaml files;
* :mod:`mpi4py`, for MPI (Message Passing Interface) functionality in Python;
* :mod:`tqdm`, for progress bars;
* :mod:`numba`, for just-in-time compilation of some Python functions;
* :mod:`palettable`, for nice colourmaps;
* :mod:`matplotlib`, for basic plotting;
* :mod:`plotly`, for advanced plotting;
* :mod:`nodejs`, for interactivity in some advanced plots;
* :mod:`ipywidgets`, for interactive plotting;
* :mod:`jupyterlab`, for convenient use of the jupyter notebooks;
* :mod:`plons`, for importing `Phantom <https://phantomsph.bitbucket.io/>`_ sph models;

All of these packages can also be found in the `conda environment file <https://github.com/Magritte-code/Magritte/blob/stable/dependencies/conda_env.yml>`_.

.. hint::

    The simplest way to setup the required python packages is using the
    `anaconda <https://www.anaconda.com/products/individual>`_ package manager.
    The Magritte conda environment can be created from the environment
    file :literal:`conda_env.yml` located in the :literal:`dependencies` directory, with

    .. code-block:: shell

        conda env create -f conda_env.yml

    This will download and install all required python packages in a newly created
    :literal:`magritte` conda environment. The environment can be activated with

    .. code-block:: shell

        conda activate magritte

    Please ensure that this environment is active whenever Magritte is compiled or used.

.. warning::

    Magritte uses plotly for some interactive plots. Plotly requires additional
    extensions to be able to render plots in a jupyter notebook or in jupyter lab. Please
    consult their `installation notes <https://plotly.com/python/getting-started/>`_ to get
    plotly working with jupyter.


Compilation
***********

Once all dependencies are in place, Magritte can be compiled.

.. hint::

    There is a shortcut script to build Magritte in the default configuration.
    From within the Magritte root directory, run:

    .. code-block:: shell

        bash build.sh

    This will create a :literal:`bin` directory in the Magritte root directory
    containing the library binary files and the executables for the tests. It will
    also create a shared object file :literal:`core.so` in the magritte python package,
    located in the :literal:`magritte` directory.

See :ref:`advanced compilation <link-advanced_compilation>` for further options.



.. _link-advanced_compilation:



Advanced compilation
********************

Compilers
=========

Corrently only the GNU gcc compiler is fully supported.
We are currently further investigating Clang and Intel compiler (:literal:`icc`) support.


GPU acceleration
================

A GPU-enabled port of Magritte to python using pytorch can be found on `GitHub <https://github.com/Magritte-code/Magritte-torch>`_.
Unless GPU acceleration is required, the C++ version of Magritte should be used, as the compiled C++ code is faster on CPU than the python version.
Not all features of the C++ version are available in the python version.
