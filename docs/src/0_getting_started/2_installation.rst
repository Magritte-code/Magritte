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

.. note::

    You may choose to not use a conda environment to install all the python dependencies for Magritte, and instead replace it with
    a regular python environment. Create a new python environment in the directory of your choice:

    .. code-block:: shell

        python -m venv magritte_env /your/environment/directory

    Then activate it (make sure it is always activated when using Magritte):

    .. code-block:: shell

        source magritte_env/bin/activate
    
    The required packages can be then be installed in your new environment using pip and the dependencies list for Magritte:

    .. code-block:: shell

        pip install -r Magritte/dependencies/requirements.txt



.. warning::

    Magritte uses plotly for some interactive plots. Plotly requires additional
    extensions to be able to render plots in a jupyter notebook or in jupyter lab. Please
    consult their `installation notes <https://plotly.com/python/getting-started/>`_ to get
    plotly working with jupyter.


Compilation
***********

Once all dependencies are in place, Magritte can be compiled. The compilation with MacOS requires extra care, it is detailed in the section :ref:`compilation on MacOS<link-macos_compilation>`.

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

.. _link-macos_compilation:

Compilation on MacOS
********************

By default, MacOS uses Clang and does not have the GNU compiler (gcc) installed. 
We do not recommend using Clang to compile Magritte because of compatibility issues with OpenMP.
Additionally, even when gcc is installed, the gcc command may still point to Clang, so a couple of 
extra steps are required to ensure that gcc is used when compiling Magritte.

If you have never installed gcc on your machine, it can be done through Homebrew 
(a package manager for MacOS, see `their webpage <https://brew.sh/>`_ for details on how to install it).

Once homebrew is installed, run the following command to install gcc:

.. code-block:: shell

    brew install gcc

gcc should now be installed, but the default gcc command may still point to Clang.
To check where the gcc command points, run the following command:

.. code-block:: shell

    gcc --version

If the output shows that the version is Apple Clang, you need manually to set the gcc command to point to the GNU compiler.
A way to do this is to set the following environment variables:

.. code-block:: shell

    export CC=/path/to/gcc/gcc-<version> 

    export CXX=/path/to/gcc/g++-<version>

Here, you should use the path to your own gcc binaries, and :literal:`<version>` is the version of gcc installed on your machine (e.g. :literal:`gcc-14` and :literal:`g++-14` for gcc.14.x.x).  
If you recompile the code often, you may want to add these two commands to your .zprofile to make them permanent.

.. hint::

    If you installed gcc through brew, its version can be found by running the following command:

    .. code-block:: shell

        brew info gcc

    cropping the version to the first number (e.g. :literal:`gcc-14` and :literal:`g++-14` for gcc.14.x.x), use the following command to find the path to your compiler:

    .. code-block:: shell

        which gcc-<version>

        which g++-<version>

    This should give you the path to your gcc binaries, which can be used in the :literal:`export` commands above. 
    The paths should look like similar to this: 
    - :literal:`/opt/homebrew/Cellar/gcc/14.2.0_1/bin/gcc-14`  
    - :literal:`/opt/homebrew/Cellar/gcc/14.2.0_1/bin/g++-14`

Once the exports are done, you can compile Magritte as described before, using:

.. code-block:: shell

    bash build.sh

We also recommend using Homebrew to install the dependencies needed by Magritte.
You can install CMake and miniconda and MPI librairies (open-mpi or mpich, and mpi4py), which are all required to compile Magritte:

.. code-block:: shell

    brew install CMake miniconda open-mpi mpi4py

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
