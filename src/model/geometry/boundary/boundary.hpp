#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


enum BoundaryCondition {Zero, Thermal, CMB};


struct Boundary
{
    Parameters parameters;

    Vector<Size> boundary2point;
    Vector<Size> point2boundary;

    Vector<BoundaryCondition> boundary_condition;
    Vector<Real>              boundary_temperature;

    BoundaryCondition set_boundary_condition (const Size b, const BoundaryCondition cd);
    BoundaryCondition get_boundary_condition (const Size b) const;

    void read  (const Io& io);
    void write (const Io& io) const;
};
