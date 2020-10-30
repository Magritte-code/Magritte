#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"


class Solver
{
    public:
        pc::multi_threading::ThreadPrivate<Vector<double>> dZ_;      ///< distance increments along the ray
        pc::multi_threading::ThreadPrivate<Vector<Size>>   nr_;      ///< corresponding point number on the ray
        pc::multi_threading::ThreadPrivate<Vector<double>> shift_;   ///< Doppler shift along the ray

        pc::multi_threading::ThreadPrivate<Vector<Real>> eta_c_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> eta_n_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> chi_c_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> chi_n_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> inverse_chi_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> tau_;

        pc::multi_threading::ThreadPrivate<Size> first_;
        pc::multi_threading::ThreadPrivate<Size> last_;
        pc::multi_threading::ThreadPrivate<Size> n_tot_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> Su_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> Sv_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> A_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> C_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> inverse_A_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> inverse_C_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> FF_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> FI_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> GG_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> GI_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> GP_;

        pc::multi_threading::ThreadPrivate<Vector<Real>> L_diag_;
        pc::multi_threading::ThreadPrivate<Matrix<Real>> L_upper_;
        pc::multi_threading::ThreadPrivate<Matrix<Real>> L_lower_;


        Size nblocks  = 512;
        Size nthreads = 512;

        Solver (const Size l, const Size w, const Size n_o_d);

        void trace (Model& model);

    // private:
        const Size length;
        const Size centre;
        const Size width;

        const Size n_off_diag;


        // void initialize (const Size l, const Size w);

        template <Frame frame>
        accel inline Size trace_ray (
            const Geometry& geometry,
            const Size      o,
            const Size      r,
            const double    dshift_max,
            const int       increment,
                  Size      id1,
                  Size      id2 );

        accel inline void set_data (
            const Size   crt,
            const Size   nxt,
            const double shift_crt,
            const double shift_nxt,
            const double dZ_loc,
            const double dshift_max,
            const int    increment,
                  Size&  id1,
                  Size&  id2 );

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

        accel inline void solve_0th_order_short_charateristics (Model& model);
        accel inline void solve_0th_order_short_charateristics (
                  Model& model,
            const Size   o,
            const Size   r,
            const double dshift_max );

        accel inline void solve_2nd_order_Feautrier (Model& model);
        accel inline void solve_2nd_order_Feautrier (
                  Model& model,
            const Size   o,
            const Size   rr,
            const Size   ar,
            const Size   f  );

        accel inline void update_Lambda (
                  Model &model,
            const Size   rr,
            const Size   f  );
};


#include "solver.tpp"
