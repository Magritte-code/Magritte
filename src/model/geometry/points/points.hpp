#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Points
{
    Parameters parameters;

    Vector <Vector3D> position;          ///< position vectors of each point
    Vector <Vector3D> velocity;          ///< velocity vectors of each point

    Vector <Size>     cum_n_neighbors;   ///< cumulative number of neighbors
    Vector <Size>         n_neighbors;   ///< number of neighbors
    Vector <Size>           neighbors;   ///< neighbors of each point

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
