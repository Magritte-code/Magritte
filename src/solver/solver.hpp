#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"


class Solver
{
    public:
        Vector<Real> dZ;      ///< distance increments along the ray
        Vector<Size> nr;      ///< corresponding point number on the ray
        Vector<Real> shift;   ///< Doppler shift along the ray

        Size nblocks  = 512;
        Size nthreads = 512;

        Solver (const Size l, const Size w);

        void trace (Model& model);
        void solve (Model& model);

    private:
        const Size length;
        const Size centre;
        const Size width;

        Vector<Size> first;
        Vector<Size> last;

        template <Frame frame>
        accel inline Size trace_ray (
            const Geometry& geometry,
            const Size      o,
            const Size      r,
            const Real      dshift_max,
            const int       increment );

        accel inline void set_data (
            const Size  crt,
            const Size  nxt,
            const Real  shift_crt,
            const Real  shift_nxt,
            const Real  dZ_loc,
            const Real  dshift_max,
            const int   increment,
                  Size& id );
};


#include "solver.tpp"