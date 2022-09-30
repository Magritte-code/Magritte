#pragma once


#include "model/model.hpp"
#include "tools/types.hpp"

///  Approximation used in the solver
/////////////////////////////////////
enum ApproximationType {None, OneLine};


struct Solver
{
    pc::multi_threading::ThreadPrivate<Vector<double>> dZ_;      ///< distance increments along the ray
    pc::multi_threading::ThreadPrivate<Vector<Size>>   nr_;      ///< corresponding point number on the ray
    pc::multi_threading::ThreadPrivate<Vector<double>> shift_;   ///< Doppler shift along the ray

    pc::multi_threading::ThreadPrivate<Vector<Real>> eta_c_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> eta_n_;

    pc::multi_threading::ThreadPrivate<Vector<Real>> chi_c_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> chi_n_;

    pc::multi_threading::ThreadPrivate<Vector<Real>> source_c_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> source_n_;

    pc::multi_threading::ThreadPrivate<Vector<Real>> inverse_chi_;

    pc::multi_threading::ThreadPrivate<Vector<Real>> tau_;

    pc::multi_threading::ThreadPrivate<Size> first_;
    pc::multi_threading::ThreadPrivate<Size> last_;
    pc::multi_threading::ThreadPrivate<Size> n_tot_;

    //indices of the sorted lines (for determining which lines to take into account when computing the opacity)
    pc::multi_threading::ThreadPrivate<Size> left_line_bound_;
    pc::multi_threading::ThreadPrivate<Size> right_line_bound_;

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

    pc::multi_threading::ThreadPrivate<Real> optical_depth_;


    Vector<Real> eta;
    Vector<Real> chi;

    Size nblocks  = 512;
    Size nthreads = 512;

    Size length;
    Size centre;
    Size width;

    Size n_off_diag;


    template <Frame frame>
    void setup (Model& model);

    void setup (const Size l, const Size w, const Size n_o_d);

    accel inline Real get_dshift_max (const Model& model, const Size o);

    template <Frame frame>
    inline void get_ray_lengths     (Model& model);

    template <Frame frame>
    inline Size get_ray_lengths_max (Model& model);


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

    template <ApproximationType approx>
    accel inline void get_eta_and_chi (
        const Model& model,
        const Size   p,
        const Size   l,
        const Real   freq,
              Real&  eta,
              Real&  chi );

    accel inline void set_sorted_line_bounds(const Model& model, Real left_freq, Real right_freq, Size p);

    template <ApproximationType approx>
    inline void compute_S_dtau_line_integrated (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext);

    inline Real compute_dtau_single_line(Model& model, Size curridx, Size nextidx, Size lineidx, Real curr_freq, Real next_freq, Real dz);

    template<ApproximationType approx>
    accel inline void compute_source_dtau (Model& model, Size currpoint, Size nextpoint, Size line, Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ, bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr, Real& Snext);

    template<ApproximationType approx>
    accel inline void update_Lambda (
              Model &model,
        const Size   rr,
        const Size   f  );


    accel inline void solve_shortchar_order_0 (Model& model);
    accel inline void solve_shortchar_order_0 (
              Model& model,
        const Size   o,
        const Size   r);
    // accel inline void solve_shortchar_order_0_ray_forward (
    //           Model& model,
    //           const Size   o,
    //           const Size   r);
    // accel inline void solve_shortchar_order_0_ray_backward (
    //           Model& model,
    //           const Size   o,
    //           const Size   r);

    // Solvers for images
    /////////////////////
    accel inline void image_feautrier_order_2 (Model& model, const Size rr);
    accel inline void image_feautrier_order_2 (Model& model, const Size o, const Size f);

    accel inline void image_feautrier_order_2_for_point     (Model& model, const Size rr, const Size p);
    accel inline void image_feautrier_order_2_for_point_loc (Model& model, const Size o,  const Size f);

    accel inline void image_optical_depth (Model& model, const Size rr);
    accel inline void image_optical_depth (Model& model, const Size o, const Size f);


    // Solvers only computing u
    ///////////////////////////
    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2 (Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2_sparse (Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2_anis (Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2 (Model& model, const Size o, const Size f);

    // Solvers for both u and v
    ///////////////////////////
    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2_uv (Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2_uv (Model& model, const Size o, const Size f);


    // Getters for emissivities, opacities, and boundary conditions
    ///////////////////////////////////////////////////////////////
    accel inline void set_eta_and_chi        (Model& model, const Size rr);
    accel inline void set_boundary_condition (Model& model) const;


    // Solvers for column densities
    ///////////////////////////////
    accel inline void set_column (Model& model) const;
    accel inline Real get_column (const Model& model, const Size o, const Size r) const;
};


#include "solver.tpp"
