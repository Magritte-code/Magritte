.. _link-installation:

Installation
############

.. note::

    This is the comprehensive installation guide. For a quick intro, see our
    :ref:`quickstart <link-quickstart>` guide.

The best way to get Magritte up-and-running is to compile it from the source code,
which can be cloned or downloaded from `GitHub <https://github.com/Magritte-code/Magritte>`_.
There are, however, a few dependencies that have to be in place before Magritte
can be compiled.
Moreover, to make Magritte useful, some python packages have to be installed which
is best done using a `conda <https://www.anaconda.com/products/individual>`_ environment.

Dependencies
************

Magritte has several dependencies, some of which are optional.

**Required**

* `CMake <https://cmake.org/>`_, version :literal:`3.18.0` or later, for building the library;

* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, version :literal:`3.3.7` or later, for some of the linear algebra;

* `pybind11 <https://github.com/pybind/pybind11>`_, version :literal:`2.2.4` or later, for binding C++ to python;

* `googletest <https://github.com/google/googletest>`_ version :literal:`1.10.0` or later, for some of the tests;

* `Paracabs <https://github.com/Magritte-code/Paracabs>`_, our custom parallelization and acceleration abstractions.

**Optional**

* `Gmsh <https://gmsh.info/>`_, version :literal:`4.6.0` or later, for meshing model geometries.



Please note that :literal:`Paracabs` might have further dependencies depending
on which paralellization and acceleration libraries are used. See
:ref:`advanced compilation <link-advanced_compilation>` for further details.


.. hint::

    There is a shortcut script to obtain the required
    dependencies: :literal:`get_dependencies.sh`. From within the :literal:`dependencies`
    directory, run:

    .. code-block:: shell

        bash get_dependencies.sh

    This will download the required dependecies, put them in the
    :literal:`dependencies` directory, and ensure the required directory sttucture
    is satisfied.

We are aware that our current dependecy management is not ideal. Any suggestions
to imporve it are welcome `here <https://github.com/Magritte-code/Magritte/issues/11>`_.


Python packages & Environment
*****************************

Magritte uses several python packages for io and to create model files.

* :mod:`numpy`, to bind the Magritte data;
* :mod:`h5py`, to read and write HDF5 data files;
* :mod:`scipy`,
* :mod:`healpy`, to sample directions from a discretized unit sphere;
* :mod:`astropy`,

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


GPU acceleration
================

Compilers
=========
