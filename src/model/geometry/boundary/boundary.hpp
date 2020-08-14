#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


struct Boundary
{
    public:
        size_t* boundary2point;
        size_t* point2boundary;
        bool*   is_on_boundary;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline void   set_npoints (const size_t n);
        inline size_t get_npoints () const;

        inline void   set_nboundary (const size_t n);
        inline size_t get_nboundary () const;

    private:
        size_t npoints;
        size_t nboundary;
};


#include "boundary.tpp"