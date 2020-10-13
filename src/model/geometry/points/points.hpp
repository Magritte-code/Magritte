#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/geometry/points/neighbors/neighbors.hpp"

<<<<<<< HEAD
//const Size nnbs = 12;
=======
// const Size nnbs = 12;
>>>>>>> 640dce3f9d1e67e0ccc99e6629edebe5685c21d5

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

<<<<<<< HEAD
//    Vector <Size> nbs; //@Frederik: what does this do?

    void read  (const Io& io);
    void write (const Io& io) const;
    //note to self: add set_neighbors and get_neighbors (should get/set all the different neighbors variables)
    //note to self: points 'deleted' should have no neighbors
=======
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
>>>>>>> 640dce3f9d1e67e0ccc99e6629edebe5685c21d5
};
