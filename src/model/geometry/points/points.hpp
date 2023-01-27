#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include <tuple>


struct Points
{
    std::shared_ptr<Parameters> parameters;

    Vector <Vector3D> position;          ///< position vectors of each point
    Vector <Vector3D> velocity;          ///< velocity vectors of each point

    Vector <Size>     cum_n_neighbors;   ///< cumulative number of neighbors
    Vector <Size>         n_neighbors;   ///< number of neighbors
    Vector <Size>           neighbors;   ///< neighbors of each point


    Points (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;

    std::tuple<Size, Size> get_neighbors(Size pointid);
};

#include "points.tpp"
