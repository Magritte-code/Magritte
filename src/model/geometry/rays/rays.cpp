#include "rays.hpp"

const string prefix = "geometry/rays/";


// Note: I can't seem to get constexpr to work correctly (linker error), thus manual template
// specialization instead Note2: It might be that pybind11 expects the manual specialization, so all
// other template functions in this file are manually specialized
template <> Size Rays ::get_direction_index<true>(const Size pointidx, const Size rayidx) const {
    return pointidx * parameters->nrays() + rayidx;
}

template <> Size Rays ::get_direction_index<false>(const Size pointidx, const Size rayidx) const {
    return rayidx;
}

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

        healpix_start_indices.resize(parameters->nrays() * parameters->npoints());
        healpix_end_indices.resize(parameters->nrays() * parameters->npoints());
        io.read_list(prefix + "healpix_start_indices", healpix_start_indices);
        io.read_list(prefix + "healpix_end_indices", healpix_end_indices);

        healpix_start_indices_sorted.resize(parameters->nrays() * parameters->npoints());
        healpix_end_indices_sorted.resize(parameters->nrays() * parameters->npoints());
        healpix_sorted_corresponding_index.resize(parameters->nrays() * parameters->npoints());

        for (Size r = 0; r < parameters->nrays() * parameters->npoints(); r++) {
            healpix_start_indices_sorted[r] = healpix_start_indices[r];
            healpix_end_indices_sorted[r] = healpix_end_indices[r];
            healpix_sorted_corresponding_index[r] = r;
        }

        // Sort for each point
        for (Size p = 0; p < parameters->npoints(); p++) {
            Size first_ray_index = get_direction_index<true>(p, 0);
            Size last_ray_index = get_direction_index<true>(p, parameters->nrays());
            std::sort(healpix_sorted_corresponding_index.dat + first_ray_index,
                      healpix_sorted_corresponding_index.dat + last_ray_index,
                      [&](Size a, Size b) {
                          return healpix_start_indices[a] < healpix_start_indices[b];
                      });
            std::sort(healpix_start_indices_sorted.dat + first_ray_index,
                      healpix_start_indices_sorted.dat + last_ray_index);
            std::sort(healpix_end_indices_sorted.dat + first_ray_index,
                      healpix_end_indices_sorted.dat + last_ray_index);
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
        io.write_list(prefix + "healpix_start_indices", healpix_start_indices);
        io.write_list(prefix + "healpix_end_indices", healpix_end_indices);
    }

    io.write_list(prefix + "weight", weight);
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

/// Returns whether the origin ray can be mapped to the target position and the resulting ray index
/// @param[in] origin_position_index: Index of the origin position
/// @param[in] origin_ray_index: Index of the origin ray direction
/// @param[in] target_position_index: Index of the target position
template <>
std::tuple<bool, Size> Rays ::get_corresponding_direction_index<false>(const Size origin_position_index, const Size origin_ray_index, const Size target_position_index) const
{
    return std::make_tuple(true, origin_ray_index);
}

// /// Returns whether the origin ray can be mapped to the target position and the resulting ray index
// /// @param[in] origin_position_index: Index of the origin position
// /// @param[in] origin_ray_index: Index of the origin ray direction
// /// @param[in] target_position_index: Index of the target position
// /// @return std::tuple(bool, Size): Whether the origin ray can be mapped to the target position and the resulting ray index \in [0, parameters.nrays()[ (which is bogus in case the bool is false)
// template <>
// std::tuple<bool, Size> Rays ::get_corresponding_direction_index<true>(const Size origin_position_index, const Size origin_ray_index, const Size target_position_index) const
// {
//     // We are comparing the healpix indices of the origin and target positions
//     Size first_ray_index = get_direction_index<true>(target_position_index, 0);
//     Size last_ray_index = get_direction_index<true>(target_position_index, parameters->nrays());
//     Size origin_ray_direction_index = get_direction_index<true>(origin_position_index, origin_ray_index);
//     Size origin_ray_start_healpix_index = healpix_start_indices[origin_ray_direction_index];
//     Size origin_ray_end_healpix_index = healpix_end_indices[origin_ray_direction_index];

//     // For the lower bound, we want the index of the first element that is less than or equal to the origin_ray_index,
//     // which is equal to the index of the first element strictly greater than the origin_ray_index minus 1.
//     auto lower_index_sorted = std::upper_bound(healpix_start_indices_sorted.dat+first_ray_index, healpix_start_indices_sorted.dat+last_ray_index, origin_ray_start_healpix_index) - 1;
//     // For the upper index, we have constructed (similar to C++ array.last()) the healpix_end_indices such that it is the index of the first element strictly greater than the origin_ray_index.
//     auto upper_index_sorted = std::lower_bound(healpix_end_indices_sorted.dat+first_ray_index, healpix_end_indices_sorted.dat+last_ray_index, origin_ray_end_healpix_index);
//     bool equal_indices = (lower_index_sorted - healpix_start_indices_sorted.dat) == (upper_index_sorted - healpix_end_indices_sorted.dat);

//     return std::make_tuple(equal_indices, healpix_sorted_corresponding_index[lower_index_sorted - healpix_start_indices_sorted.dat] - first_ray_index);
// }

/// Returns whether the origin ray can be mapped to the target position and the resulting ray index
/// @param[in] origin_position_index: Index of the origin position
/// @param[in] origin_ray_index: Index of the origin ray direction
/// @param[in] target_position_index: Index of the target position
/// @return std::tuple(bool, Size): Whether the origin ray can be mapped to the target position and the resulting ray index \in [0, parameters.nrays()[ (which is bogus in case the bool is false)
template <>
std::tuple<bool, Size> Rays ::get_corresponding_direction_index<true>(const Size origin_position_index, const Size origin_ray_index, const Size target_position_index) const
{
    // We are comparing the healpix indices of the origin and target positions
    Size first_ray_index = get_direction_index<true>(target_position_index, 0);
    Size last_ray_index = get_direction_index<true>(target_position_index, parameters->nrays());
    Size origin_ray_direction_index = get_direction_index<true>(origin_position_index, origin_ray_index);
    Size origin_ray_start_healpix_index = healpix_start_indices[origin_ray_direction_index];
    Size origin_ray_end_healpix_index = healpix_end_indices[origin_ray_direction_index];

    // For the lower bound, we want the index of the first element that is less than or equal to the origin_ray_index,
    // which is equal to the index of the first element strictly greater than the origin_ray_index minus 1.
    auto lower_index_sorted = std::upper_bound(healpix_start_indices_sorted.dat+first_ray_index, healpix_start_indices_sorted.dat+last_ray_index, origin_ray_start_healpix_index) - 1;
    // For the upper index, we have constructed (similar to C++ array.last()) the healpix_end_indices such that it is the index of the first element strictly greater than the origin_ray_index.
    auto upper_index_sorted = std::lower_bound(healpix_end_indices_sorted.dat+first_ray_index, healpix_end_indices_sorted.dat+last_ray_index, origin_ray_end_healpix_index);
    bool equal_indices = (lower_index_sorted - healpix_start_indices_sorted.dat) == (upper_index_sorted - healpix_end_indices_sorted.dat);

    //DEBUG: test if limiting to up to 1 layer up works; assumes that the weights of the finer layers are equal to 1/4 of the coarser layer
    Real weight_origin_ray = get_weight<true>(origin_position_index, origin_ray_index);
    const Size target_ray_index = healpix_sorted_corresponding_index[lower_index_sorted - healpix_start_indices_sorted.dat] - first_ray_index;
    Real weight_target_ray = get_weight<true>(target_position_index, healpix_sorted_corresponding_index[lower_index_sorted - healpix_start_indices_sorted.dat] - first_ray_index);
    // std::cout<<"origin_ray_index: "<<origin_ray_index<<" target ray index: "<<target_ray_index<<std::endl;
    // std::cout<<"origin_ray_weight: "<<weight_origin_ray<<" target ray weight: "<<weight_target_ray<<std::endl;
    if (weight_origin_ray < 1.0/2.0 * weight_target_ray) {
    // if (weight_origin_ray < 1.0/5.0 * weight_target_ray) {
        // std::cout<<"false"<<std::endl;
        equal_indices = false;
    }
    // equal_indices = false;

    // return std::make_tuple(equal_indices, healpix_sorted_corresponding_index[lower_index_sorted - healpix_start_indices_sorted.dat] - first_ray_index);
    return std::make_tuple(equal_indices, target_ray_index);
}