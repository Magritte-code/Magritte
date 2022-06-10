#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Turbulence
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model

    Vector<Real> vturb2;                      ///< [.] microturbulence over c all squared


    Turbulence (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;
};
