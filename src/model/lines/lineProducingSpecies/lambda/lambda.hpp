#pragma once


#include "tools/types.hpp"


class Lambda
{
    public:
        Real3 Ls;    ///< values
        Size3 nr;    ///< position indices

        Real1 Lss;   ///< linearized values
        Size1 nrs;   ///< linearized position indices

        Size1 size;


        inline void initialize (const Parameters &parameters, const size_t nrad_new);

        inline void clear ();

        inline void linearize_data ();

        inline int MPI_gather ();


        inline Size index_first (const Size p, const Size k) const;
        inline Size index_last  (const Size p, const Size k) const;

        inline Real get_Ls (const Size p, const Size k, const Size index) const;
        inline Size get_nr (const Size p, const Size k, const Size index) const;

        inline Size get_size (const Size p, const Size k) const;

        inline void add_element (const Size p, const Size k, const Size nr, const Real Ls);

    private:
        Size npoints;
        Size nrad;
};


#include "lambda.tpp"
