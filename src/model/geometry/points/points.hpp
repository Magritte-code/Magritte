#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/geometry/points/neighbors/neighbors.hpp"



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

    // Vector <Size> nbs;

    void read  (const Io& io);
    void write (const Io& io) const;

    void print()
    {
        for (Size r = 0; r < 10; r++)
        {
            position[r].print();
        }
    }
};
