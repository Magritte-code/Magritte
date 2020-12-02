Getting started
###############

.. Warning::
    The instructions below are currently only tested on Linux.


Instalation
***********

Dependencies
============


There is a shortcut script to install the required dependecies. From within the
Magritte root directory run:

.. code-block:: shell

    bash dependencies/get_dependencies.sh


Compilation
===========

Once all dependencies are in place, the default version of Magritte can be compiled with:

.. code-block:: shell

    bash build.sh

This will create a :literal:`bin` directory in the Magritte root directory
containing the library binary files and the executables for the tests. It will
also create a shared object file :literal:`core.so` in the magritte python package,
located in the :literal:`magritte` directory.


Usage
*****

The `examples <https://github.com/Magritte-code/Examples>`_ repository contains some jupyter notebooks and equivalent python
scripts working examples
Here we only highlight some peculiarities.
Currently Magritte is not yet a proper python package.


Examples
========

Some working examples can be found in `this`
