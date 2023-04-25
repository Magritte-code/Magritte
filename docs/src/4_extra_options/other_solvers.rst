Other solvers
=============

Magritte currently contains several different solvers for computing NLTE radiative transfer. In this document, we will explain their usecases.

The default solver recommended by the examples is the 2nd order Feautrier solver (Feautrier 1964).
As it is both accurate and fast, it is a fine option to default to.

.. code-block:: python

   magritte.core.compute_level_populations(use_Ng_acceleration, max_niterations)

The memory usage of this solver can be too high for larger models. By no longer storing the intensities for each direction,
the memory cost is reduced significantly. The sparse solver can be called similarly using:

.. code-block:: python

   magritte.core.compute_level_populations_sparse(use_Ng_acceleration, max_niterations)

Magritte also offers the option to compute the intensities using the formal solution.
This solver can be used to compute the directional intensities.

.. code-block:: python

   magritte.core.compute_level_populations_shortchar(use_Ng_acceleration, max_niterations)
