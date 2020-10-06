#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/geometry/points/neighbors/neighbors.hpp"

//const Size nnbs = 12;

struct Points
{
    Parameters parameters;

    Vector <Vector3D> position;          ///< position vectors of each point
    Vector <Vector3D> velocity;          ///< velocity vectors of each point

    Neighbors curr_neighbors;   ///Does everything the old data 'structure' did
//    Vector <Size>     cum_n_neighbors;   ///< cumulative number of neighbors
//    Vector <Size>         n_neighbors;   ///< number of neighbors each point has
//    Vector <Size>           neighbors;   ///< neighbors of each point, listed after eachother
    ///< e.g.: [neighbors of point 1, neighbors of point 2, ...]

//    Vector <Size> nbs; //@Frederik: what does this do?

    void read  (const Io& io);
    void write (const Io& io) const;
    //note to self: add set_neighbors and get_neighbors (should get/set all the different neighbors variables)
    //note to self: points 'deleted' should have no neighbors
};
