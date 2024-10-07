NLTE in very low density regimes
================================

In low density regimes, the level populations might start to maser. 
This means that the resulting line opacity (and optical depth) can become close to zero or even negative.
Unfortunately, our solvers cannot handle a rapid change in source function (emissivity/opacity) too well, and might produce unphysical results.
Therefore we both implement a minimum value for the line opacity and an interpolation scheme to avoid discontinuities in the line opacity.

This is why we implement a workaround by setting the minimum allowed value for the line opacity (starting from version 0.7.0).
This value is set to 1e-13 [W sr−1 m−3] by default, but can be changed by the user. Higher values can overestimate the impact of the less dense model regions on the intensity. Lower values can impact the number of interpolation points required for the computations (see below).

.. code-block:: python

    parameters.min_line_opacity = 1e-13#default

In versions 0.6.0 until 0.9.0, the value was set to 1e-10 [W sr−1 m−3] by default, and the interpolation procedure was not yet implemented.

----------------------------------------------------------------------------------------------------------------------------

In NLTE line radiative transfer, the line opacity can change orders of magnitude between successive positions, leading to numerical artefacts.
We implement an interpolation scheme to avoid these discontinuities in the line opacity.
To regulate the number of interpolation points, the user can set the following parameter:

.. code-block:: python

    parameters.max_interpolation_diff = 1.4#default

This parameter ensures that the (normalized) line opacities at max differ a factor *max_interpolation_diff* between two successive points.