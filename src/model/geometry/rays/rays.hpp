#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Rays
{
    Parameters parameters;

    Vector<Vector3D> direction;
    Vector<Size>     antipod;
    Vector<Real>     weight;

    void read  (const Io& io);
    void write (const Io& io) const;
};


#include "rays.tpp"