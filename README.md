<img src="docs/src/images/Magritte_logo_plain.svg" alt="logo" width="350"/>

[![Build status](https://github.com/Magritte-code/Magritte/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/Magritte-code/Magritte/actions/workflows/build-and-test.yml)
[![Documentation Status](https://readthedocs.org/projects/magritte/badge/?version=stable)](https://magritte.readthedocs.io/en/stable/?badge=stable)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03905/status.svg)](https://doi.org/10.21105/joss.03905)
---

Welcome to the Magritte repository! Magritte is a open-source software library for
simulating radiation transport, developed at
[University College London](https://www.ucl.ac.uk/) (UCL, UK) and
[KU Leuven](https://www.kuleuven.be/english/) (Belgium).

Magritte is currently mainly used for post-processing hydrodynamical simulations of
astrophysical models by creating synthetic observations, but the techniques could
in principle also be applied more general (e.g. plasma and nuclear physics).  

It can either be used as a Python package or as a C++ library.
Magritte uses a deterministic ray-tracer with a formal solver that currently focusses on
line radiative transfer (see
[De Ceuster et al. 2019](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D/abstract)
for more details). The plot below shows synthetic observations made with
Magritte at three different inclinations for a hydro simulation of the stellar wind around
an asymptotic giant branch (AGB) star as it is perturbed by a companion.

<img src="docs/src/_static/movie.gif" alt="movie"/>

The hydro model was created by Jan Bolte using [MPI-AMRVAC](http://amrvac.org/). See the
[examples](https://magritte.readthedocs.io/en/latest/1_examples/index.html) in the
documentation to learn how these synthetic observations were created with Magritte.


## Documentation
Please find our online documentation [here](https://magritte.readthedocs.io).
In particular, see our
[quickstart documentation](https://magritte.readthedocs.io/en/latest/0_getting_started/0_quickstart.html)
for a quick intro.


## Papers about Magritte
The following list of papers might provide further insights in the inner workings of
Magritte:
* _Magritte II: Adaptive ray-tracing, mesh construction and reduction_
([arXiv](https://arxiv.org/abs/2011.14998), [MNRAS](https://doi.org/10.1093/mnras/staa3199));
* _Magritte I: Non-LTE atomic and molecular line modelling_
([arXiv](https://arxiv.org/abs/1912.08445),
[MNRAS](https://doi.org/10.1093/mnras/stz3557));
* _3D Line Radiative Transfer & Synthetic Observations with
Magritte_ ([JOSS](https://doi.org/10.21105/joss.03905)).

Please note that some features presented in these papers might not yet be implemented
and documented in the latest release of Magritte.


## Issues & Contact
Please report any issues with Magritte or its documentation
[here](https://github.com/UCL/Magritte/issues).
If you need any further help, please contact
[Frederik De Ceuster](https://freddeceuster.github.io/).


## Cite
Please contact the authors of the papers referenced above if you want to use
Magritte in your research. We are currently working on documentation and
examples to facilitate its independent use. Until then, please
[contact us](https://freddeceuster.github.io/).


## Developers & Contributors
**Developers**
* Frederik De Ceuster
* Thomas Ceulemans

**Scientific & Technical advisors**
* Ward Homan
* Jan Bolte
* Jeremy Yates
* Leen Decin
* Peter Boyle
* James Hetherington

**Contributors**
* Silke Maes
* Jolien Malfait
* Atulit Srivastava
* Mats Esseldeurs
* Arnout Coenegrachts

## Acknowledgements
FDC was for most of the development supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC, JB, WH, and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.
TC is a PhD fellow of the Research Foundation - Flanders (FWO).
