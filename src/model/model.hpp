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

        inline void set_npoints (const Size npoints);
        inline Size get_npoints () const;

    private:
        Size npoints;
};


#include "model.tpp"