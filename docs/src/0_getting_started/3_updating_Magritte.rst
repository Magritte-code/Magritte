.. _link-updating_Magritte:

Updating Magritte
#################

This part of the documentation assumes that one has a previous version of magritte installed
and the core (non-python, non-submodule) dependencies of Magritte did not change. Otherwise, check out the :ref:`prerequisites <link-prerequisites>`.

Pull
****

To get the latest magritte version, we pull from the repository, including all submodules.

.. code-block:: shell

    git pull --recurse-submodules


Update dependencies
*******************

We update the conda environment, as the python dependencies might have changed.

.. code-block:: shell

    conda activate magritte
    conda env update -f dependencies/conda_env.yml


Build
*****

Magritte can be build using its default configuration.

.. code-block:: shell

    bash build.sh


The latest version of magritte is now installed on your machine.
