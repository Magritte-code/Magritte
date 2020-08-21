#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"

const Size nnbs = 12;

class Points
{
    public:
        Vector <Vector3D> position;          ///< position vectors of each point
        Vector <Vector3D> velocity;          ///< velocity vectors of each point

        Vector <Size>     cum_n_neighbors;   ///< cumulative number of neighbors
        Vector <Size>         n_neighbors;   ///< number of neighbors
        Vector <Size>           neighbors;   ///< neighbors of each point

        Vector <Size> nbs;

//        Vector* position;          ///< position vectors of each point
//        Vector3D* velocity;          ///< velocity vectors of each point

//        size_t* cum_n_neighbors;   ///< cumulative number of neighbors
//        size_t*     n_neighbors;   ///< number of neighbors
//        size_t*       neighbors;   ///< neighbors of each point

//         Points();
//        ~Points();
//
        inline void set_npoints (const Size n);
        inline Size get_npoints () const;

        void read  (const Io& io);
        void write (const Io& io) const;

    private:
        Size npoints;
        Size tot_n_neighbors;
};


#include "points.tpp"