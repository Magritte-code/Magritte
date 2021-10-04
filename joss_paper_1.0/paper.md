---
title: '3D Line Radiative Transfer & Synthetic Observations with \textsc{Magritte}'
tags:
  - Python
  - C++
  - astronomy
  - radiative transfer
authors:
  - name: Frederik De Ceuster^[Corresponding author.]
    orcid: 0000-0001-5887-8498
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Thomas Ceulemans
    #orcid: 0000-0001-5887-8498
    affiliation: "2"
  - name: Atulit Srivastava
    #orcid: 0000-0003-3997-9642
    affiliation: "2"
  - name: Jeremy Yates
    orcid: 0000-0003-1954-8749
    affiliation: "2"
  - name: Leen Decin
    orcid: 0000-0002-5342-8612
    affiliation: "1, 3"
  - name: Peter Boyle
    orcid: 0000-0002-8960-1587
    affiliation: "4, 5"
  - name: James Hetherington
    orcid: 0000-0001-6993-0319
    affiliation: "6, 7"
affiliations:
 - name: Department of Physics & Astronomy, University College London, London, UK
   index: 1
 - name: Institute of Astronomy, KU Leuven, Leuven, Belgium
   index: 2
 - name: School of Chemistry, University of Leeds, Leeds, UK
   index: 3
 - name: School of Physics & Astronomy, The University of Edinburgh, Edinburgh, UK
   index: 4
 - name: Brookhaven National Laboratory, NY, USA
   index: 5
 - name: Department of Computer Science, University College London, London, UK
   index: 6
 - name: The Alan Turing Institute, London, UK
   index: 7
date: 1 October 2021
bibliography: paper.bib


# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary

\textsc{Magritte} is a software library for 3D radiative transfer and synthetic observations,
which currently focusses on atomic and molecular lines.
It is mainly written in C++ and can either be used as a Python package or
as a C++ library. To compute the radiation field, a deterministic ray-tracer
and a formal solver are employed, i.e. rays are traced through the model after which the
radiative transfer equation is solved along those rays `[@DeCeuster:2019]`. This
is in contrast to most radiative transfer solvers which employ (probabilistic)
Monte Carlo techniques `[@Noebauer:2019]`. By vitrue of minimal assumptions about
the underlying geometric structure of a model, Magritte can handle structured and
unstructured input meshes, as well as smoothed-particle hydrodynamics (SPH) data, and,
futhermore, provides tools to optimize different input meshes for radiative transfer
`[@DeCeuster:2020]`.


# Statement of need

Electromagnetic radiation is a key component in many astrophysical simulations.
Not only does it dictate what we can or cannot observe, it can provide e.g. radiation
pressure, efficient heating and cooling mechanisms, and opens up a range of new
chemical pathways due to photo-reactions. Magritte can be used as a
general-purpose radiative transfer solver, but was particularly designed for
line radiative transfer in complex 3D morphologies, such as, for instance,
encountered in the stellar winds around evolved stars `[see @Decin:2021]`.
Recent high-resolution observations exposed the intricate and intrinsically 3D
morphologies of these objects `[@Decin2020]`. The sheer amount of complexity that is
observed, makes it difficult to interpret the observations and necessitates the use of
3D hydrodynamical/chemical models to study their origin and evolution `[@ElMellah:2020;
@Bolte2021; @Maes2021; @Malfait2021]`. Their intricate morpho-kinematics, moreover,
makes their appearance in observations far from trivial. Therefore, it is essential to
understand how these models would appear in observations. This can be achieved, by
making synthetic observations with Magritte (see e.g. \autoref{fig:example}).

![Example of a synthetic observation of the CO($v=0$, $J=1-0$)-transition, made with Magritte for a hydrodynamics model of an asymptotic giant branch (AGB) star, as it is perturbed by a companion. This is model \textsc{v10e50} in `[@Malfait:2021]` \label{fig:example}](example.png)


# Future work

Currently Magritte is mainly being used for post-processing hydrodynamical/chemical
simulations by creating synthetic observations of the models such that they can be
compared with real observations . However, Magritte is still under active development.
In future work, we aim to clear the path for on-the-fly radiative transfer in those
simulations. Aside from being a practical research tool, Magritte also aims to be a
foundation for research in computational radiative transfer. Current active research
topics include: efficient parallelisation and acceleration strategies on modern
high-performance computing systems, acceleration of convergence in the non-linear
coupling between radiation field and medium, and uncertainty quantification
in radiative transfer through probabilistic numerical methods.


# Acknowledgements

FDC is supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.

# References
