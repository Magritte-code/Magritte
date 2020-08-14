#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"


class Solver
{
    public:
        double* dZ;      ///< distance increments along the ray
        size_t* nr;      ///< corresponding point number on the ray
        double* shift;   ///< Doppler shift along the ray

    private:
        const size_t length;
        const size_t centre;
        const size_t width;

        size_t* first;
        size_t* last;

         Solver (const size_t l, const size_t w);
        ~Solver ();


        void trace (const Model& model);

        template <Frame frame>
        inline void trace_ray (
            const Geometry& geometry,
            const size_t    o,
            const size_t    r,
            const double    dshift_max,
            const int       increment   );

        inline void set_data (
            const size_t  crt,
            const size_t  nxt,
            const double  shift_crt,
            const double  shift_nxt,
            const double  dZ_loc,
            const double  dshift_max,
            const int    increment,
                  size_t& id         );
};


#include "solver.tpp"