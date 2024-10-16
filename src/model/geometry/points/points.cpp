#include "points.hpp"

#include <assert.h>

const string prefix = "geometry/points/";

void Points ::read(const Io& io) {
    cout << "Reading points..." << endl;

    parameters->set_npoints(io.get_length(prefix + "position"));
    parameters->set_npoints(io.get_length(prefix + "velocity"));

    position.resize(parameters->npoints());
    velocity.resize(parameters->npoints());

    Double2 position_buffer(parameters->npoints(), Double1(3));
    Double2 velocity_buffer(parameters->npoints(), Double1(3));

    io.read_array(prefix + "position", position_buffer);
    io.read_array(prefix + "velocity", velocity_buffer);

    Real min_x = std::numeric_limits<Real>::max();
    Real max_x = std::numeric_limits<Real>::min();
    Real min_y = std::numeric_limits<Real>::max();
    Real max_y = std::numeric_limits<Real>::min();
    Real min_z = std::numeric_limits<Real>::max();
    Real max_z = std::numeric_limits<Real>::min();

    for (Size p = 0; p < parameters->npoints(); p++) {
        position[p] = Vector3D(position_buffer[p][0], position_buffer[p][1], position_buffer[p][2]);

        velocity[p] = Vector3D(velocity_buffer[p][0], velocity_buffer[p][1], velocity_buffer[p][2]);

        if (position[p].x() < min_x) min_x = position[p].x();
        if (position[p].x() > max_x) max_x = position[p].x();
        if (position[p].y() < min_y) min_y = position[p].y();
        if (position[p].y() > max_y) max_y = position[p].y();
        if (position[p].z() < min_z) min_z = position[p].z();
        if (position[p].z() > max_z) max_z = position[p].z();
    }

    center = Vector3D((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2);

    const Size totnnbs = io.get_length(prefix + "neighbors");

    cum_n_neighbors.resize(parameters->npoints());
    n_neighbors.resize(parameters->npoints());
    neighbors.resize(totnnbs);

    io.read_list(prefix + "n_neighbors", n_neighbors);
    io.read_list(prefix + "neighbors", neighbors);

    cum_n_neighbors[0] = 0;

    for (Size p = 1; p < parameters->npoints(); p++) {
        cum_n_neighbors[p] = cum_n_neighbors[p - 1] + n_neighbors[p - 1];
    }

    position.copy_vec_to_ptr();
    velocity.copy_vec_to_ptr();

    cum_n_neighbors.copy_vec_to_ptr();
    n_neighbors.copy_vec_to_ptr();
    neighbors.copy_vec_to_ptr();
}

void Points ::write(const Io& io) const {
    cout << "Writing points..." << endl;

    Double2 position_buffer(parameters->npoints(), Double1(3));
    Double2 velocity_buffer(parameters->npoints(), Double1(3));

    for (size_t p = 0; p < parameters->npoints(); p++) {
        position_buffer[p] = {position[p].x(), position[p].y(), position[p].z()};
        velocity_buffer[p] = {velocity[p].x(), velocity[p].y(), velocity[p].z()};
    }

    io.write_array(prefix + "position", position_buffer);
    io.write_array(prefix + "velocity", velocity_buffer);

    io.write_list(prefix + "n_neighbors", n_neighbors);
    io.write_list(prefix + "neighbors", neighbors);
}
