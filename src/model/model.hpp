#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"


class Model
{
    public:
        Geometry geometry;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline void   set_npoints (const size_t npoints);
        inline size_t get_npoints () const;

    private:
        size_t npoints;
};


#include "model.tpp"