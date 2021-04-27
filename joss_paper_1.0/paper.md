---
title: 'Magritte: a modern software library for radiative transfer and synthetic observations'
tags:
  - Python
  - C++
  - astronomy
  - radiative transfer
authors:
  - name: Frederik De Ceuster^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: 0000-0001-5887-8498
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
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
 - name: Department of Physics and Astronomy, University College London, Gower Place, London, WC1E 6BT, UK
   index: 1
 - name: Department of Physics and Astronomy, Institute of Astronomy, KU Leuven, Celestijnenlaan 200D, 3001 Leuven, Belgium
   index: 2
 - name: School of Chemistry, University of Leeds, Leeds LS2 9JT, UK
   index: 3
 - name: School of Physics and Astronomy, The University of Edinburgh, Edinburgh EH9 3FD, UK
   index: 4
 - name: Brookhaven National Laboratory, Upton, NY 11973, USA
   index: 5
 - name: Department of Computer Science, University College London, Bloomsburry, London, WC1E 6EA, UK
   index: 6
 - name: The Alan Turing Institute, 96 Euston Road, Kings Cross, London, NW1 2DB, UK
   index: 7
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Radiative processes play a key role in various throughout astronomy, physics, and engineering.
Various simulations throughout astronomy, physics and engineering require an accurate model for
the transport of radiation.


# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements

FDC is supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.

# References
