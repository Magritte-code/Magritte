#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/geometry/points/multiscale/multiscale.hpp"


struct Points
{
    Parameters parameters;

    Vector <Vector3D> position;          ///< position vectors of each point
    Vector <Vector3D> velocity;          ///< velocity vectors of each point

    Multiscale multiscale;               ///< data structure containing the neighbor data at each level

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
