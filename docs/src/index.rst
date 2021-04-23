Magritte documentation
######################

.. Warning::
    These pages are still under construction!

Welcome to the Magritte documentation! Magritte is an open-source software
library for simulating radiation transport, developed at `University College
London <https://www.ucl.ac.uk/>`_ (UCL, UK) and `KU Leuven
<https://www.kuleuven.be/english/>`_ (Belgium).

Magritte is currently mainly used for post-processing hydrodynamical simulations by
creating synthetic observations. The plot below shows synthetic observations made with
Magritte at three different inclinations for a hydro simulation of the stellar wind around
an asymptotic giant branch (AGB) star as it is perturbed by a companion.

.. raw:: html

    <video width="100%" controls playsinline autoplay muted loop>
      <source src="_static/movie.webm">
      Your browser does not support the video tag.
    </video>

The hydro model was created by Jan Bolte using `MPI-AMRVAC <http://amrvac.org/>`_. See
the :ref:`examples <link-examples>` to learn how the synthetic observations were created
with Magritte.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   0_getting_started/index
   1_examples/index
   2_benchmarks/index
   3_python_api_documentation/index
   4_cpp_api_documentation/index


Papers about Magritte
*********************

The following list of papers might provide further insights in the inner workings of Magritte:

* Magritte: **Adaptive ray-tracing, mesh construction and reduction**,
  *F. De Ceuster, J. Bolte, W. Homan, S. Maes, J. Malfait, L. Decin, J. Yates, P. Boyle, J. Hetherington*, 2020
  (`arXiv <https://arxiv.org/abs/2011.14998>`_,
  `MNRAS <https://doi.org/10.1093/mnras/staa3199>`_);


* Magritte: **Non-LTE atomic and molecular line modelling**,
  *F. De Ceuster, W. Homan, J. Yates, L. Decin, P. Boyle, J. Hetherington*, 2019
  (`arXiv <https://arxiv.org/abs/1912.08445>`_,
  `MNRAS <https://doi.org/10.1093/mnras/stz3557>`_);

Please note that some features presented in these papers might not yet be implemented
and documented in the latest release of Magritte.


Contact
*******

Please report any `issues <https://github.com/Magritte-code/Magritte/issues>`_ with
Magritte or its documentation `here <https://github.com/Magritte-code/Magritte/issues>`_.
If you need any further help, please contact `Frederik De Ceuster
<https://www.kuleuven.be/wieiswie/en/person/00101884>`_.


Developers & Contributors
*************************

**Developers**

* Frederik De Ceuster
* Thomas Ceulemans
* Atulit Srivastava

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


Acknowledgements
****************

FDC is supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC, JB, WH, and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.
