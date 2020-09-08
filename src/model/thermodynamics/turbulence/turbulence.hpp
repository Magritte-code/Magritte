#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Turbulence
{
    Parameters parameters;

    Vector<Real> vturb2;   ///< [.] microturbulence over c all squared

    void read  (const Io& io);
    void write (const Io& io) const;
};
