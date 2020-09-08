#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Boundary
{
    Parameters parameters;

    Vector<Size> boundary2point;
    Vector<Size> point2boundary;

    void read  (const Io& io);
    void write (const Io& io) const;
};


#include "boundary.tpp"