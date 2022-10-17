Approximations in Magritte
==========================

Currently, the most expensive computation in Magritte involves computing sums
over the line opacities :math:`χ_{ij}(ν)` and line emissivities :math:`η_{ij}(ν)`:

.. math::
   χ_{ij}(ν) &= hν/4π (n_jB_{ji}-n_iB_{ij})ϕ_{ij}(ν)\\
   η_{ij}(ν) &= hν/4π n_iA_{ij}ϕ_{ij}(ν)


However, very few of the profile functions :math:`ϕ_{ij}(ν)` will be nonzero at any time during the computation.
Currently, Magritte takes advantage of this by first computing which lines lie close enough to significantly contribute.
Afterwards, it will only sum over the nearby lines.

To disable this feature (not recommended; only for comparison purposes), one can call

.. code-block:: python

   parameters.sum_opacity_emissivity_over_all_lines = True #(default = False)

A further improvement can be made if all lines are spaced far enough from each other.
In that case, we do not need to search which lines lie close anymore and can instead 'sum' over a single line.
In order for this approximation to be valid, no lines should overlap (including doppler shift effects).

To enable this feature, one can call

.. code-block:: python

   parameters.one_line_approximation = True #(default = False)

.. Warning::
   Do not set both parameters.sum_opacity_emissivity_over_all_lines and parameters.one_line_approximation to True.
   As these settings are contradictory, this results in undefined behavior.
