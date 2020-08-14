#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


class Points
{
    public:
        Vector* position;          ///< position vectors of each point
        Vector* velocity;          ///< velocity vectors of each point

        size_t* cum_n_neighbors;   ///< cumulative number of neighbors
        size_t*     n_neighbors;   ///< number of neighbors
        size_t*       neighbors;   ///< neighbors of each point

//         Points();
//        ~Points();
//
        inline void   set_npoints (const size_t n);
        inline size_t get_npoints () const;

        void read  (const Io& io);
        void write (const Io& io) const;

    private:
        size_t npoints;
        size_t tot_n_neighbors;
};


#include "points.tpp"