#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

struct Rays {
    std::shared_ptr<Parameters> parameters;

    Vector<Vector3D> direction; // Direction Can either have dimensions [Nrays, 3] or
                                // [Nrays*Npoints, 3], in order to keep backward compatibility
    Vector<Size> antipod;       // Opposite direction index; assume to be the same for all points;
                                // dimensions [Nrays]
    Vector<Real> weight;        // Weight of each ray; dimensions [Nrays, 3] or [Nrays*Npoints, 3]
    bool use_adaptive_directions =
        false; // Whether to use a different set of directions for each ray

    Rays(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;

    Size get_direction_index(
        const Size pointidx, const Size rayidx) const; // linearized direction index
    Vector3D get_direction(const Size pointidx, const Size rayidx) const;
    Vector3D get_antipod(const Size pointidx, const Size rayidx) const;
    Size get_antipod_index(const Size rayidx) const;
    Real get_weight(const Size pointidx, const Size rayidx) const;
};
