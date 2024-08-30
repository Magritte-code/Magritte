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
    Vector<Size> healpix_start_indices; // Start of indices at finest directional discretization corresponding to each ray direction;
                //optional (only for adaptive directions); dimensions [Nrays*Npoints]
    Vector<Size> healpix_end_indices; // End (+1; same meaning as vector.last() in C++) of indices at finest directional discretization corresponding to each ray direction;
                //optional (only for adaptive directions); dimensions [Nrays*Npoints]
                //e.g.; we have directions corresponding to [[0], [1], [2], [3], [4,5,6,7]] (for a single point)
                //then healpix_start_indices = [0, 1, 2, 3, 4] and healpix_end_indices = [1, 2, 3, 4, 8]origin_ray_end_healpix_index
    Vector<Size> healpix_start_indices_sorted;
    Vector<Size> healpix_end_indices_sorted;
    Vector<Size> healpix_sorted_corresponding_index;//For converting the sorted indices back to the original indices

    Rays(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;

    template <bool use_adaptive_directions>
    Size get_direction_index(
        const Size pointidx, const Size rayidx) const; // linearized direction index
    template <bool use_adaptive_directions>
    Vector3D get_direction(const Size pointidx, const Size rayidx) const;
    template <bool use_adaptive_directions>
    Vector3D get_antipod(const Size pointidx, const Size rayidx) const;
    Size get_antipod_index(const Size rayidx) const;
    template <bool use_adaptive_directions>
    Real get_weight(const Size pointidx, const Size rayidx) const;
    template <bool use_adaptive_directions>
    std::tuple<bool, Size> get_correspoding_direction_index(const Size origin_position_index, const Size origin_ray_index, const Size target_position_index) const;
};
