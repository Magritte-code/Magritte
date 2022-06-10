#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Rays
{
    std::shared_ptr<Parameters> parameters;

    Vector<Vector3D> direction;
    Vector<Size>     antipod;
    Vector<Real>     weight;


    Rays (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;
};
