NLTE in very low density regimes
================================

In low density regimes, the level populations might start to maser. 
This means that the resulting line opacity (and optical depth) can become close to zero or even negative.
Unfortunately, our solvers cannot handle this situation too well, and might produce unphysical results.

In more recent versions of the code (starting from version 0.5.3), a workaround has been implemented to avoid this situation by resetting specific levels at specific places to LTE.
This might work in some cases (with the caveat of some artifacts being produced in the resulting spectrum), 
but might as well produce new artifacts in when doing NLTE line radiative transfer in extremely low density regions. 

To restore the previous behavior (without the new workaround), one can call

.. code-block:: python

    parameters.population_inversion_fraction = -1