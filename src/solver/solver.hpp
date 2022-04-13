#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"
#include "tools/interpolation.hpp"
#include <set>


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


        //comoving approach
        // pc::multi_threading::ThreadPrivate<Vector<Size>> n_rays_through_point_;//unsigned char might be a bit too short
        //technically, both next vars can be threadprivate, if we could guarantee that the same ray distribution would be used? (between initialization and actually using it)
        Vector<Vector<Size>> points_to_trace_ray_through;
        // Vector<Size> n_points_to_trace_ray_through;
        Matrix<Size> n_rays_through_point;//ray, point
        // pc::multi_threading::ThreadPrivate<Size> n_rays_to_trace_;
        //I will assume no more than 2^16-1 rays go through a specific point if choosing unsigned int
        //However, just go with a Size for absolute safety (max upper bound on number rays traced can be assumed to be the number of points)
        pc::multi_threading::ThreadPrivate<Matrix<Real>> J_;//for purpose of collecting J in all different threads
        pc::multi_threading::ThreadPrivate<Vector<unsigned char>> real_pt_;
        pc::multi_threading::ThreadPrivate<Vector<Real>> temp_intensity_;

        // pc::multi_threading::ThreadPrivate<Vector<Real>> boundary_frequencies;
        // pc::multi_threading::ThreadPrivate<Vector<Real>> boundary_frequencies;

        // Real dshift_max_lower_bound;


        // Kernel approach
        Vector<Real> eta;
        Vector<Real> chi;

        // SparseMatrix<Real> covariance;
        // Matrix<Real> L2_kernel_p;

        Size nblocks  = 512;
        Size nthreads = 512;

        // Solver () {};
        // Solver (const Size l, const Size w, const Size n_o_d);

        template <Frame frame>
        void setup (Model& model);
        void setup (const Size l, const Size w, const Size n_o_d);

        template <Frame frame>
        void static_setup(Model& model);

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

        template <Frame frame>
        accel inline Size trace_ray_static (
            const Geometry& geometry,
            const Size      o,
            const Size      r,
            const double    dshift_max,
            const int       increment,
                  Size      id1,
                  Size      id2 );


        accel inline void trace_ray_points (
            const Geometry& geometry,
            const Size      o,
            const Size      rdir,
            const Size      rsav);

        inline void get_static_rays_to_trace (Model& model);

        inline void compute_u_static (Model& model);

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

        accel inline void set_data_static (
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

        accel inline void get_eta_and_chi_static (
            const Model& model,
            const Size   p,
            const Real   freq,
            const double shift,
                  Real&  eta,
                  Real&  chi ) const;

        accel inline void update_Lambda (
                  Model &model,
            const Size   rr,
            const Size   f  );


        accel inline void solve_shortchar_order_0 (Model& model);
        accel inline void solve_shortchar_order_0 (
                  Model& model,
            const Size   o,
            const Size   r,
            const double dshift_max );

        accel inline void solve_shortchar_static (Model& model);
        accel inline void solve_shortchar_static (
                  Model& model,
            const Size   o,
            const Size   r,
            const double dshift_max );

        // accel inline void solve_shortchar_order_0_comoving (Model& model);

        accel inline void solve_feautrier_order_2 (Model& model);
        accel inline void solve_feautrier_order_2 (
                  Model& model,
            const Size   o,
            const Size   rr,
            const Size   ar,
            const Size   f  );

        accel inline void image_feautrier_order_2 (Model& model, const Size rr);
        accel inline void image_feautrier_order_2 (
                  Model& model,
            const Size   o,
            const Size   rr,
            const Size   ar,
            const Size   f  );

        accel inline void image_feautrier_order_2_for_point (Model& model, const Size rr, const Size p);
        accel inline void image_feautrier_order_2_for_point (
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
