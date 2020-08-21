#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


class Rays
{
    public:
        Vector<Vector3D> direction;
        Vector<Size>     antipod;
        Vector<Real>     weight;
//        Vector3D* direction;
//        size_t*   antipod;
//        double*   weight;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline Size get_nrays () const;

    private:
        Size nrays;
};


#include "rays.tpp"