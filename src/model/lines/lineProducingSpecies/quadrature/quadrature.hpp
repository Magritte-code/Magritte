#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


struct Quadrature
{
    public:
        Double1 roots;
        Double1 weights;

        void read  (const Io &io, const Size l);
        void write (const Io &io, const Size l) const;

        accel inline void set_nquads (const Size n);
        accel inline Size get_nquads () const;

    private:
        size_t nquads;   ///< number frequency quadrature points
};


#include "quadrature.tpp"