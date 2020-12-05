Usage
#####

The `examples repository <https://github.com/Magritte-code/Examples>`_  contains
jupyter notebooks (and equivalent python scripts) for research examples using
Magritte. On the :ref:`examples page <link-examples>` we showcase some of the
results. Here we only highlight some basics to get you started.


Currently, Magritte is not yet a proper python package. Hence, every python script
using Magritte should add it to the :literal:`PYTHONPATH` with:

.. code-block:: python

    from sys import path
    path.append('path/to/magritte/root/directory')
