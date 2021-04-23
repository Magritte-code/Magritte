.. _link-post-processing:


Post-processing hydro snapshots
###############################

We show how Magritte can be used to post-process snapshots of hydrodynamics simulations.
We consider snapshots of both `MPI-AMRVAC <http://amrvac.org/>`_, a solver based on
adaptive mesh refinement (AMR, `Xia et al. 2018
<https://ui.adsabs.harvard.edu/abs/2018ApJS..234...30X/abstract>`_), and `Phantom
<https://phantomsph.bitbucket.io/>`_ which uses smoothed-particle hydrodynamics (SPH,
`Price et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018PASA...35...31P/abstract>`_).

.. Note::
    For general examples on how to generate Magritte models from analytic or tabulated
    data, please refer to the :ref:`creating model examples <link-creating_models>`.

First we show how to build Magritte models for these snapshots and then we demonstrate
how to reduce them to smaller approximate models for fast radiative transfer
(`De Ceuster et al. 2020
<https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5194D/abstract>`_). 


.. toctree::
   :maxdepth: 1
   :caption: AMRVAC snapshots:
   
   0_create_AMRVAC_1D.ipynb
   1_create_AMRVAC_3D.ipynb
   2_reduce_AMRVAC_3D.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Phantom snapshots:  

   3_create_Phantom_3D.ipynb
   4_reduce_Phantom_3D.ipynb

