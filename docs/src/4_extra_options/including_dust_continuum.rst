Including dust or continuum sources in Magritte
===============================================

Starting from Magritte version 0.9.11, you can include additional continuum sources in your model.
These are treated in LTE; therefore you need to provide a dust/continuum opacity table (per point, per frequency) and temperature (per point).
This can be done using:

.. code-block:: python

    model.dust.dust_opacities.set(continuum_opacity_table)# in 1/m dims: [N_points, N_continuum_frequencies]
    model.dust.dust_frequencies.set(frequency_table)# in Hz dims: [N_continuum_frequencies]
    model.dust.dust_temperatures.set(temperature_table)# in K dims: [N_points]

Evidently, we have also provides some tools to help calculate the dust opacities from tabulated data.
In particular, if you have a table from the Jena dust database (https://www2.astro.uni-jena.de/Laboratory/OCDB/), you can use the following commands:

.. code-block:: python

    import magritte.tools as tools
    continuum_reference_density = astropy.Quantity # in mass density [e.g. g/cm^3]; given by the continuum data; used for rescaling the values
    frequency_grid, complex_refractive_index = tools.read_dust_opacity_table(continuum_file_location, ASTROPY_UNIT_OF_FREQ_WAVELENGTH_...)
    # reads the file and returns the frequency grid (first column), multiplied whatever astropy unit you specified, and the complex refractive index (unitless, third column) as a numpy array
    frequency_grid_Hz, absorption_coeff_per_density = tools.convert_dust_opacity_table_to_SI(frequency_grid, complex_refractive_index, continuum_ref_density)
    # returns the frequency grid in Hz and absorption coeff in m^2/kg (needs to be multiplied with the mass density to get the absorption coeff per m)
    #for this example, assume we scale the continuum/dust number density with the H2 number density
    continuum_opacity_table = (absorption_coeff_per_density/1000) * fraction_density_dust * nH2 * 2.016 / 6.022e23 #H2: 2.016 g/mol, 6.022e23 particles/mol -> particles/m3 to g/m3

In which the last step converts the absorption coefficient per mass density (m^2/kg) to an absorption coefficient per number density (1/m).
In this release, we also included options to image the model using a custom frequency grid (given that equidistant grids might not be optimal for imaging the continuum).

.. code-block:: python

    # Instead of model.compute_spectral_discretization(...)
    model.set_custom_spectral_discretization(frequency_grid_Hz)# np.array, in Hz, dims: [N_frequencies_to_image]
    model.compute_image_new(...) #as usual


.. Warning::

    Even though you might have a model for which you only have continuum sources (no line sources), you still need to include a dummy species producing lines.
    You can even set its abundance to zero, but Magritte needs at least one species to properly initialize a model. See tests/benchmarks/analytic/continuum_implementation.py for an example.