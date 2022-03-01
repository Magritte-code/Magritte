---
title: '3D Line Radiative Transfer & Synthetic Observations with Magritte'
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
  - name: Ward Homan
    orcid: 0000-0001-7314-5081
    affiliation: "3"
  - name: Jan Bolte
    orcid: 0000-0002-7991-9663
    affiliation: "2"
  - name: Jeremy Yates
    orcid: 0000-0003-1954-8749
    affiliation: "1"
  - name: Leen Decin
    orcid: 0000-0002-5342-8612
    affiliation: "2, 4"
  - name: Peter Boyle
    orcid: 0000-0002-8960-1587
    affiliation: "5, 6"
  - name: James Hetherington
    orcid: 0000-0001-6993-0319
    affiliation: "7"
affiliations:
 - name: Department of Physics and Astronomy, University College London, London, UK
   index: 1
 - name: Institute of Astronomy, KU Leuven, Leuven, Belgium
   index: 2
 - name: Institut d’Astronomie et d’Astrophysique, Université Libre de Bruxelles, Brussels, Belgium
   index: 3
 - name: School of Chemistry, University of Leeds, Leeds, UK
   index: 4
 - name: School of Physics and Astronomy, The University of Edinburgh, Edinburgh, UK
   index: 5
 - name: Brookhaven National Laboratory, NY, USA
   index: 6
 - name: Department of Computer Science, University College London, London, UK
   index: 7
date: 1 October 2021
bibliography: paper.bib


# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary

Electromagnetic radiation is a key component in many astrophysical simulations.
Not only does it dictate what we can or cannot observe, it can provide radiation
pressure, efficient heating and cooling mechanisms, and opens up a range of new
chemical pathways due to photo-reactions. Magritte is a software library that can be
used as a general-purpose radiative transfer solver, but was particularly designed for
line radiative transfer in complex 3D morphologies, such as, for instance,
encountered in the stellar winds around evolved stars [see @Decin:2021].
It is mainly written in C++ and can either be used as a Python package or
as a C++ library. To compute the radiation field, a deterministic ray-tracer
and a formal solver are employed, i.e., rays are traced through the model and the
radiative transfer equation is solved along those rays [@DeCeuster:2019]. This
is in contrast to most radiative transfer solvers which employ (probabilistic)
Monte Carlo techniques [@Noebauer:2019]. By virtue of minimal assumptions about
the underlying geometric structure of a model, Magritte can handle structured and
unstructured input meshes, as well as smoothed-particle hydrodynamics (SPH) data.
Furthermore, tools are provided to optimise different input meshes for radiative transfer
[@DeCeuster:2020].


# Statement of need

Recent high-resolution observations exposed the intricate and intrinsically 3D
morphologies of stellar winds around evolved stars [@Decin:2020]. The sheer amount of complexity that is
observed, makes it difficult to interpret the observations and necessitates the use of
3D hydrodynamics, chemistry and radiative transfer models to study their origin and
evolution [@ElMellah:2020; @Maes:2021; @Malfait:2021]. Their intricate
morpho-kinematics, moreover, makes their appearance in (synthetic) observations far from evident
(see e.g.\ the intricate structures in \autoref{fig:example}). Therefore, to study these and other complex
morpho-kinematical objects, it is essential to understand how their models would
appear in observations. This can be achieved, by creating synthetic observations
with Magritte.
Examples and analytic as well as cross-code benchmarks can be found in the documentation and in [@DeCeuster:2019; @DeCeuster:2020].

![Example of a synthetic observation of the CO($v=0$, $J=1-0$) transition, created with Magritte for a hydrodynamics model of an asymptotic giant branch (AGB) star, as it is perturbed by a companion [this is model \textsc{v10e50} in @Malfait:2021]. \label{fig:example}](example.png)


# Future work

Currently, Magritte is mainly being used for post-processing hydrodynamics
simulations by creating synthetic observations, such that the models can be
compared with real observations. However, Magritte is still under active development.
In future work, we aim to explicitly include continuum radiation and improve on the computational speed to clear the path for on-the-fly radiative transfer in hydrodynamics models. Furthermore, aside from being a
practical research tool, we also aim for Magritte to be the starting point for
further research in computational radiative transfer. Current active research topics
include: efficient parallelisation and acceleration strategies on modern
high-performance computing systems, acceleration of convergence in the non-linear
coupling between the radiation field and the medium, and uncertainty quantification in
radiative transfer through probabilistic numerical methods [e.g., @DeCeuster:2021].


# Acknowledgements

FDC is supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC, JB, and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.
TC is a PhD fellow of the Research Foundation – Flanders (FWO).
WH acknowledges support from the Fonds de la Recherche Scientifique (FNRS) through grant 40000307.

# References
