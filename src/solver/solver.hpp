#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"


class Solver
{
    public:
        Vector<double> dZ;      ///< distance increments along the ray
        Vector<Size>   nr;      ///< corresponding point number on the ray
        Vector<double> shift;   ///< Doppler shift along the ray



        pc::multi_threading::ThreadPrivate<Vector<Real>> eta_crt;
        pc::multi_threading::ThreadPrivate<Vector<Real>> eta_nxt;

        pc::multi_threading::ThreadPrivate<Vector<Real>> chi_crt;
        pc::multi_threading::ThreadPrivate<Vector<Real>> chi_nxt;

        pc::multi_threading::ThreadPrivate<Vector<Real>> drho;
        pc::multi_threading::ThreadPrivate<Vector<Real>> dtau;

        pc::multi_threading::ThreadPrivate<Vector<Real>> tau;


        Size nblocks  = 512;
        Size nthreads = 512;

        Solver (const Size l, const Size w);

        void trace (Model& model);
        void solve (Model& model);

    // private:
        const Size length;
        const Size centre;
        const Size width;

        Vector<Size> first;
        Vector<Size> last;

        // void initialize (const Size l, const Size w);

        template <Frame frame>
        accel inline Size trace_ray (
            const Geometry& geometry,
            const Size      o,
            const Size      r,
            const double    dshift_max,
            const int       increment );

        accel inline void set_data (
            const Size   crt,
            const Size   nxt,
            const double shift_crt,
            const double shift_nxt,
            const double dZ_loc,
            const double dshift_max,
            const int    increment,
                  Size&  id );

        accel inline Real gaussian (const Real width, const Real diff) const;
        accel inline Real planck   (const Real temp,  const Real freq) const;

        accel inline Real boundary_intensity (
            const Model& model,
            const Size   p,
            const Real   freq ) const;

        accel inline void get_eta_and_chi (
            const Model& model,
            const Size   p,
            const Real   freq,
                  Real&  eta,
                  Real&  chi );

        accel inline void solve_0th_order_short_charateristics (
                  Model& model,
            const Size   o,
            const Size   r,
            const double dshift_max );
};


#include "solver.tpp"
