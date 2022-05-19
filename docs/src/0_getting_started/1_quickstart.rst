.. _link-quickstart:

Quickstart
##########

.. Warning::
    Please ensure that the :ref:`prerequisites <link-prerequisites>` are in place before continuing with the installation.

.. Warning::
    This is the (too) quick intallation guide. For details, see our
    :ref:`comprehensive installation <link-installation>` guide.


Clone
*****

Magritte has to be compiled from its source code, which can be cloned using:

.. code-block:: shell

    git clone --recursive https://github.com/Magritte-code/Magritte.git

from our `GitHub <https://github.com/Magritte-code/Magritte>`_ repository.


Setup
*****

The Magritte `conda <https://www.anaconda.com/products/individual>`_ environment
can be created from an environment file with:

.. code-block:: shell

    conda env create -f dependencies/conda_env.yml

which installs all required packages in the :literal:`magritte` conda
environment, and can be activated with:

.. code-block:: shell

    conda activate magritte

This environment has to be active whenever Magritte is compiled or used!

.. Note::
    If you didn't configure your shell with the Anaconda installation, the above command won't work, but the same result can be achieved with:

    .. code-block:: shell

        source activate magritte


Build
*****

Magritte can easily be build in its default configuration with:

.. code-block:: shell

    bash build.sh

.. Hint:: By default, the c- and c++-compilers to use for building, get determined by

    .. code-block:: shell
    
        which gcc
        which g++

    One can optionally define the environment variables CC, and CXX to point respectively to a c- and c++-compiler.
    This is mostly useful for MacOs users, as gcc by default points to Clang, which has some issues with OpenMP.

This will create the necessary files for the Magritte python package.


Run
***

If all the above worked, go download and experiment with some of our :ref:`examples
<link-examples>`!
