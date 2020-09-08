#pragma once


#include "io/io.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
//    Chemistry      chemistry;
//    Thermodynamics thermodynamics;

    void read  (const Io& io);
    void write (const Io& io) const;
};


#include "model.tpp"