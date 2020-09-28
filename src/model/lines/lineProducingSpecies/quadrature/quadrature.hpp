#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Quadrature
{
    Parameters parameters;

    Vector<Real> roots;
    Vector<Real> weights;

    void read  (const Io& io, const Size l);
    void write (const Io& io, const Size l) const;
};
