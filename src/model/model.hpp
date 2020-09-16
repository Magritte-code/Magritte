#pragma once


#include "io/io.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;

    void read  (const Io& io);
    void write (const Io& io) const;
};


#include "model.tpp"