#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Quadrature
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    Vector<Real> roots;                       ///< Quadrature roots
    Vector<Real> weights;                     ///< Quadrature weights


    Quadrature (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io, const Size l);
    void write (const Io& io, const Size l) const;
};
