.. _link-quickstart:

Quickstart
##########

.. Warning::
    This is the (too) quick intallation guide. For details, see our
    :ref:`comprehensive installation <link-installation>` guide.


Get Magritte
************

Magritte has to be compiled from its source code, which can be either cloned

.. code-block:: shell

    git clone https://github.com/Magritte-code/Magritte.git

or downloaded

.. code-block:: shell

    wget https://github.com/Magritte-code/Magritte/archive/master.zip
    unzip master.zip
    mv Magritte-master Magritte

from `GitHub <https://github.com/Magritte-code/Magritte>`_. The
:literal:`Magritte` directory this gives you will be the Magritte root directory.


Get dependencies
****************

The required dependencies can be obtained with

.. code-block:: shell

    bash dependencies/get_dependencies.sh

This will download the required dependecies and put them in the
:literal:`dependencies` directory.


Setup python environment
************************

The Magritte `conda <https://www.anaconda.com/products/individual>`_ environment can be created from the environment file

.. code-block:: shell

    conda env create -f dependencies/conda_env.yml

which installs all required packages in the :literal:`magritte` conda
environment, and can be activated with

.. code-block:: shell

    conda activate magritte

This environment has to be active whenever Magritte is compiled or used!


Build Magritte
**************

Magritte can be build in the default configuration with

.. code-block:: shell

    bash build.sh

This will create the necessary files for the Magritte python package.


Run Magritte
************

If all the above worked, go download and experiment with the
`examples <https://github.com/Magritte-code/Examples>`_!
