#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

enum BoundaryCondition { Zero, Thermal, CMB };

struct Boundary {
    std::shared_ptr<Parameters> parameters;

    Vector<Size> boundary2point; // maps boundary index ∈[0, nboundary-1] to point index ∈[0, npoints-1]
    Vector<Size> point2boundary; // maps point index to boundary index. Maps non-boundary points to
                                 // parameters->npoints() (invalid value)

    Vector<BoundaryCondition> boundary_condition;
    Vector<Real> boundary_temperature;

    Boundary(std::shared_ptr<Parameters> params) : parameters(params){};

    BoundaryCondition set_boundary_condition(const Size b, const BoundaryCondition cd);
    BoundaryCondition get_boundary_condition(const Size b) const;

    void read(const Io& io);
    void write(const Io& io) const;
};
