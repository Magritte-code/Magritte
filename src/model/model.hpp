#pragma once


#include "io/io.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"
#include "radiation/radiation.hpp"


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;
    Radiation      radiation;

    enum SpectralDiscretisation {None, SD_Lines, SD_Image}
         spectralDiscretisation = None;

    void read  (const Io& io);
    void write (const Io& io) const;

    int compute_inverse_line_widths     ();
    int compute_spectral_discretisation ();
    int compute_spectral_discretisation (const Real width);
};
