#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


class Rays
{
    public:
        Vector* direction;
        size_t* antipod;
        double* weight;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline size_t get_nrays () const;

    private:
        size_t nrays;
};


#include "rays.tpp"