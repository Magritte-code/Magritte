#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


struct Boundary
{
    public:
        Vector<Size> boundary2point;
        Vector<Size> point2boundary;
//        Vector<bool>   is_on_boundary;

//        size_t* boundary2point;
//        size_t* point2boundary;
//        bool*   is_on_boundary;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline void set_npoints (const Size n);
        inline Size get_npoints () const;

        inline void set_nboundary (const Size n);
        inline Size get_nboundary () const;

    private:
        Size npoints;
        Size nboundary;
};


#include "boundary.tpp"