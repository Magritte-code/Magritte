.. _link-quickstart:

Quickstart
##########

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

The Magritte `conda <https://www.anaconda.com/products/individual>`_ environment can be created from an environment file with:

.. code-block:: shell

    conda env create -f dependencies/conda_env.yml

which installs all required packages in the :literal:`magritte` conda
environment, and can be activated with:

.. code-block:: shell

    conda activate magritte

This environment has to be active whenever Magritte is compiled or used!


Build
*****

Magritte can easily be build in its default configuration with:

.. code-block:: shell

    bash build.sh

This will create the necessary files for the Magritte python package.


Run
***

If all the above worked, go download and experiment with some of our
`examples <https://github.com/Magritte-code/Examples>`_!
