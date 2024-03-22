.. _link-prerequisites:

Prerequisites
#############

Throughout this documentation, we assume a Unix-based operating system (e.g. Linux or MacOS).
For Windows users, we recommend to use the Windows Subsystem Linux (WSL; see e.g. `here <https://www.windowscentral.com/install-windows-subsystem-linux-windows-10>`_).

Below is a list of the software required to be able to compile and run Magritte.

* `GCC <https://gcc.gnu.org/>`_ (version :literal:`5.0.0` or later): to compile the C++ part of Magritte.
* `CMake <https://cmake.org/>`_ (version :literal:`3.18.0` or later): for building Magritte, organising compilation and linking;
* `Anaconda <https://www.anaconda.com/blog/individual-edition-2020-11>`_ (optional): for managing the required Python packages;
* `Gmsh <https://gmsh.info/>`_ (optional): for constructing geometric meshes for the models;

Although Anaconda is optional and you could, in priciple, use any other way of installing the required Python dependencies, we recommend and assume the use of Anaconda.
The list of required Python packages can be found in the `conda environment file <https://github.com/Magritte-code/Magritte/blob/stable/dependencies/conda_env.yml>`_.
Also Gmsh is optional, however, it will be used in most examples in this documentation to create models.

Once these prerequisites are in place, Magritte can be installed following the :ref:`quickstart <link-quickstart>`.

.. Warning::
    On MacOS, gcc by default points to Apple Clang (which is NOT the GNU compiler). Apple Clang won't work, due to an issue with OpenMP.
    Hence, GNU gcc has to be installed separately, e.g. using `Homebrew <https://brew.sh/>`_ or `Macports <https://www.macports.org/>`_.


.. Hint::
    Although we recommend to first try to install the required software in your preferred usual way, the folder `Magritte/dependencies/ <https://github.com/Magritte-code/Magritte/tree/stable/dependencies>`_ contains a set of bash scripts: :literal:`get_conda.sh`, :literal:`get_cmake.sh`,
    and :literal:`get_gmsh.sh`, to quickly obtain the appropriate versions of miniconda3 (for Anaconda), CMake, and
    Gmsh respectively. However, these scripts don't properly install the software, but rather download an executable binary, which is useful, e.g. if you don't have root access.
    
    To use these, first clone the Magritte `GitHub <https://github.com/Magritte-code/Magritte>`_ repository
    
    .. code-block:: shell

        git clone --recursive https://github.com/Magritte-code/Magritte.git
    
    enter the newly created Magritte directory, and from there run the relevant script.
    
    .. code-block:: shell
    
        bash dependencies/get_conda.sh
        
    .. code-block:: shell
    
        bash dependencies/get_cmake.sh
        
    .. code-block:: shell
    
        bash dependencies/get_gmsh.sh

    ⚠️  Don't forget to include the relevant paths to your shell's :literal:`$PATH` variable!
    (For convenience, the scripts will print the relevant path.)

.. Note::
    We are working on a set of pre-compiled versions of Magritte that would greatly simply the installation process.
    `Let us know <https://github.com/Magritte-code/Magritte/issues>`_, if you have any suggestions or want to help!
