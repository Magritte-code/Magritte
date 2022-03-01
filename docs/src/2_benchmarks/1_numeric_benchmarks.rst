Numeric benchmarks
##################


Van Zadelhoff et al. (2002)
***************************

The `van Zadelhoff et al. (2002) <https://ui.adsabs.harvard.edu/abs/2002A%26A...395..373V>`_
benchmark comprises two main benchmark problems for non-LTE line radiative transfer.
The results for Magritte on this benchmark and others are discussed and compared with
`RATRAN <https://personal.sron.nl/~vdtak/ratran/frames.html>`_
and `LIME <https://github.com/lime-rt/lime>`_ in
`De Ceuster et al. (2019) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D>`_.
The scripts to setup and run these benchmarks are provided with Magritte and can be found `here <https://github.com/Magritte-code/Magritte/tree/stable/tests/benchmarks/numeric>`_.


Problem 1
=========

The first problem considers a fictitious two-level species in a spherically symmetric
cloud, without velocity field, with a constant temperature distribution, and a
quadratucally decaying density distribution. A line data file for the fictitious
two-level species (in `LAMDA <https://home.strw.leidenuniv.nl/~moldata/>`_ format)
is provided with Magritte and can be found
`here <https://github.com/Magritte-code/Magritte/blob/master/tests/data/test.txt>`_.


Problem 2
=========

The second problem has a more realistic setup and considers the lines of HCO+ in a snapshot of an inside-out collapse model.
The model consists of 50 logarithmically spaced grid points. In each grid point the radial velocity, gas temperature, micro-turbulence, and both HCO+ and H2 abundances are given.
This data was provided on the benchmark  website (which is offline now), but can still be found `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/benchmarks/numeric/vanZadelhoff_2a.in>`_ and `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/benchmarks/numeric/vanZadelhoff_2b.in>`_.
The line data file for HCO+ comes from the LAMDA data base and can be found `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/data/hco%2B.txt>`_