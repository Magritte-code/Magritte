Getting started
###############

.. Warning::
    The instructions below were only tested under Linux.


Installation
************

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

The `examples repository <https://github.com/Magritte-code/Examples>`_  contains
jupyter notebooks (and equivalent python scripts) for research examples using
Magritte. On the :ref:`examples page <link-examples>` we showcase some of the
results. Here we only highlight some basics to get you started.


Currently, Magritte is not yet a proper python package. Hence, every python script
using Magritte should add it to the :literal:`PYTHONPATH` with:

.. code-block:: python

    from sys import path
    path.append('path/to/magritte/root/directory')
