Usage
#####

The `examples repository <https://github.com/Magritte-code/Examples>`_  contains
jupyter notebooks (and equivalent python scripts) for research examples using
Magritte. On the :ref:`examples page <link-examples>` we showcase some of the
results. Here we only highlight some basics to get you started.

As already stated in the quickstart and installation guide, we highly recommend
to use Magritte within a `conda <https://www.anaconda.com/products/individual>`_
environment. The Magritte conda environment, in which all examples should work,
can be created from an environment file with:

.. code-block:: shell

    conda env create -f dependencies/conda_env.yml

which installs all required packages, and can be activated with:

.. code-block:: shell

    conda activate magritte

Please note that when using this conda environment, Magritte has to be installed
in it, i.e. this environment has to be active whenever Magritte is compiled!
