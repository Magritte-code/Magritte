#include "rays.hpp"

const string prefix = "geometry/rays/";

void Rays ::read(const Io& io) {
    cout << "Reading rays..." << endl;

    Size len_dir =
        io.get_length(prefix + "direction"); // length of first axis of the direction array

    if (len_dir == parameters->nrays()) {
        parameters->set_hnrays(parameters->nrays() / 2);
        use_adaptive_directions = false;
    } else if (len_dir == parameters->nrays() * parameters->npoints()) {
        parameters->set_hnrays(parameters->nrays() / 2);
        use_adaptive_directions = true;
    } else {
        parameters->set_nrays(
            len_dir); // will error, and give some more info on the sizes of the arrays
    }

    antipod.resize(parameters->nrays());

    if (!use_adaptive_directions) {
        weight.resize(parameters->nrays());
        io.read_list(prefix + "weight", weight);

        direction.resize(parameters->nrays());
        Double2 direction_buffer(parameters->nrays(), Double1(3));

        io.read_array(prefix + "direction", direction_buffer);

        for (Size r = 0; r < parameters->nrays(); r++) {
            direction[r] =
                Vector3D(direction_buffer[r][0], direction_buffer[r][1], direction_buffer[r][2]);
        }
    } else { // use adaptive directions == true
        weight.resize(parameters->nrays() * parameters->npoints());
        io.read_list(prefix + "weight", weight);

        direction.resize(parameters->nrays() * parameters->npoints());
        Double2 direction_buffer(parameters->nrays() * parameters->npoints(), Double1(3));

        io.read_array(prefix + "direction", direction_buffer);

        for (Size r = 0; r < parameters->nrays() * parameters->npoints(); r++) {
            direction[r] =
                Vector3D(direction_buffer[r][0], direction_buffer[r][1], direction_buffer[r][2]);
        }
    }

    const double tolerance = 1.0E-9;

    // We assume the antipod index to remain the same for all points
    for (Size r1 = 0; r1 < parameters->nrays(); r1++) {
        for (Size r2 = 0; r2 < parameters->nrays(); r2++) {
            if ((direction[r1] + direction[r2]).squaredNorm() < tolerance) {
                antipod[r1] = r2;
            }
        }
    }

    direction.copy_vec_to_ptr();
    antipod.copy_vec_to_ptr();
    weight.copy_vec_to_ptr();
}

void Rays ::write(const Io& io) const {
    cout << "Writing rays..." << endl;

    if (!use_adaptive_directions) {
        Double2 direction_buffer(parameters->nrays(), Double1(3));

        for (Size r = 0; r < parameters->nrays(); r++) {
            direction_buffer[r] = {direction[r].x(), direction[r].y(), direction[r].z()};
        }
        io.write_array(prefix + "direction", direction_buffer);
    } else {
        Double2 direction_buffer(parameters->nrays() * parameters->npoints(), Double1(3));

        for (Size r = 0; r < parameters->nrays() * parameters->npoints(); r++) {
            direction_buffer[r] = {direction[r].x(), direction[r].y(), direction[r].z()};
        }
        io.write_array(prefix + "direction", direction_buffer);
    }

    io.write_list(prefix + "weight", weight);
}

// Note: I can't seem to get constexpr to work correctly (linker error), thus manual template
// specialization instead Note2: It might be that pybind11 expects the manual specialization, so all
// other template functions in this file are manually specialized

template <> Size Rays ::get_direction_index<true>(const Size pointidx, const Size rayidx) const {
    return pointidx * parameters->nrays() + rayidx;
}

template <> Size Rays ::get_direction_index<false>(const Size pointidx, const Size rayidx) const {
    return rayidx;
}

/// Get the direction of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Vector3D Rays ::get_direction<true>(const Size pointidx, const Size rayidx) const {
    return direction[get_direction_index<true>(pointidx, rayidx)];
}

/// Get the direction of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Vector3D Rays ::get_direction<false>(const Size pointidx, const Size rayidx) const {
    return direction[get_direction_index<false>(pointidx, rayidx)];
}

/// Get the antipodal direction of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Vector3D Rays ::get_antipod<true>(const Size pointidx, const Size rayidx) const {
    return direction[get_direction_index<true>(pointidx, get_antipod_index(rayidx))];
}

/// Get the antipodal direction of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Vector3D Rays ::get_antipod<false>(const Size pointidx, const Size rayidx) const {
    return direction[get_direction_index<false>(pointidx, get_antipod_index(rayidx))];
}

/// Get the antipodal direction index
/// @param[in] rayidx: Index of the ray
Size Rays ::get_antipod_index(const Size rayidx) const { return antipod[rayidx]; }

/// Get the weight of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Real Rays ::get_weight<true>(const Size pointidx, const Size rayidx) const {
    return weight[get_direction_index<true>(pointidx, rayidx)];
}

/// Get the weight of a ray
/// @param[in] pointidx: Index of the point
/// @param[in] rayidx: Index of the ray
template <> Real Rays ::get_weight<false>(const Size pointidx, const Size rayidx) const {
    return weight[get_direction_index<false>(pointidx, rayidx)];
}
