#include "geometry.hpp"

void Geometry ::read(const Io& io) {
    points.read(io);
    rays.read(io);
    boundary.read(io);

    lengths.resize(parameters->hnrays(), parameters->npoints());
}

void Geometry ::write(const Io& io) const {
    points.write(io);
    rays.write(io);
    boundary.write(io);
}
