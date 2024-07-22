NLTE in very low density regimes
================================

In low density regimes, the level populations might start to maser. 
This means that the resulting line opacity (and optical depth) can become close to zero or even negative.
Unfortunately, our solvers cannot handle a rapid change in source function (emissivity/opacity) too well, and might produce unphysical results.

This is why we implement a workaround by setting the minimum allowed value for the line opacity (starting from version 0.7.0).
This value is set to 1e-10 [W sr−1 m−3] by default, but can be changed by the user. Higher values might overestimate the impact of the less dense model regions on the intensity, while lower value can lead to numerical artifacts in the NLTE intensities.

.. code-block:: python

    parameters.min_line_opacity = 1e-10#default

In versions 0.5.3 until 0.6.0, the code implemented another workaround to avoid this situation by resetting specific levels at specific places to LTE.
This might work in some cases (with the caveat of some artifacts being produced in the resulting spectrum), 
but might as well produce new artifacts in when doing NLTE line radiative transfer in extremely low density regions. 

To restore that behavior (without the new workaround), one can call

.. code-block:: python

    parameters.population_inversion_fraction = 1.01#previous default
    parameters.min_line_opacity = -1

Finally, to restore the behavior of the code before version 0.5.3 (not recommended when using NLTE), one can call

.. code-block:: python

    parameters.population_inversion_fraction = -1
    parameters.min_line_opacity = -1