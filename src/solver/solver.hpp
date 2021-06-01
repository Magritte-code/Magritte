#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"


class Solver
{
    public:
        VectorTP<double> dZ;      ///< distance increments along the ray
        VectorTP<Size>   nr;      ///< corresponding point number on the ray
        VectorTP<double> shift;   ///< Doppler shift along the ray

        VectorTP<Real> eta_c;
        VectorTP<Real> eta_n;

        VectorTP<Real> chi_c;
        VectorTP<Real> chi_n;

        VectorTP<Real> inverse_chi;

        VectorTP<Real> tau;

        TP<Size> first;
        TP<Size> last;
        TP<Size> n_tot;

        VectorTP<Real> Su;
        VectorTP<Real> Sv;

        VectorTP<Real> A;
        VectorTP<Real> C;
        VectorTP<Real> inverse_A;
        VectorTP<Real> inverse_C;

        VectorTP<Real> FF;
        VectorTP<Real> FI;
        VectorTP<Real> GG;
        VectorTP<Real> GI;
        VectorTP<Real> GP;

        VectorTP<Real> L_diag;
        MatrixTP<Real> L_upper;
        MatrixTP<Real> L_lower;


        // Kernel approach
        Vector<Real> eta;
        Vector<Real> chi;

        // SparseMatrix<Real> covariance;
        // Matrix<Real> L2_kernel_p;

        Size nblocks  = 1; //512;
        Size nthreads = 1; //512;

        // Solver () {};
        // Solver (const Size l, const Size w, const Size n_o_d);

        template <Frame frame>
        void setup (Model& model);
        void setup (const Size l, const Size w, const Size n_o_d);

        accel inline Real get_dshift_max (
            const Model& model,
            const Size   o     );

        template <Frame frame>
        inline void get_ray_lengths     (Model& model);
        template <Frame frame>
        inline Size get_ray_lengths_max (Model& model);


    // private:
        Size length;
        Size centre;
        Size width;

        Size n_off_diag;


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
                  Real&  chi ) const;

        accel inline void update_Lambda (
                  Model &model,
            const Size   rr,
            const Size   f  );


              inline void solve_shortchar_order_0 (Model& model);
        accel inline void solve_shortchar_order_0 (
                  Model& model,
            const Size   o,
            const Size   r,
            const double dshift_max );

              inline void solve_feautrier_order_2 (Model& model);
        accel inline void solve_feautrier_order_2 (
                  Model& model,
            const Size   o,
            const Size   rr,
            const Size   ar,
            const Size   f  );

              inline void image_feautrier_order_2 (Model& model, const Size rr);
        accel inline void image_feautrier_order_2 (
                  Model& model,
            const Size   o,
            const Size   rr,
            const Size   ar,
            const Size   f  );


        accel inline Real     kernel (const Vector3D d) const;
        accel inline Real     kernel (const Model& model, const Size r, const Size p1, const Size p2) const;
        accel inline Real  L1_kernel (const Model& model, const Size r, const Size p1, const Size p2) const;
        accel inline Real  L2_kernel (const Model& model, const Size r, const Size p1, const Size p2) const;
        accel inline Real L12_kernel (const Model& model, const Size r, const Size p1, const Size p2) const;

        accel inline void solve_kernel_method (
                  Model& model,
            const Size   r,
            const Size   f );

        accel inline void set_eta_and_chi        (Model& model) const;
        accel inline void set_boundary_condition (Model& model) const;
};


#include "solver.tpp"
