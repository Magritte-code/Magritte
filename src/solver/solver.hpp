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

    pc::multi_threading::ThreadPrivate<Real> optical_depth_;

    // Comoving approach


    pc::multi_threading::ThreadPrivate<Vector<Size>> freqsplits_;//resize to nlines as overestimate for size of this
    pc::multi_threading::ThreadPrivate<Size> n_freqsplits_;



    Vector<Vector<Size>> points_to_trace_ray_through;
    Matrix<Size> n_rays_through_point;//raydir, point ; might not be used that much in the future
    //hmm, we need to weight the rays somehow; for now we just only use the closest ray to determine intensity
    Matrix<Size>> closest_ray;//raydir, point (contains closest rayid?)
    Matrix<Size>> min_ray_distsqr;//raydir, pointid
    // pc::multi_threading::ThreadPrivate<Size> n_rays_to_trace_;
    //I will assume no more than 2^16-1 rays go through a specific point if choosing unsigned int
    //However, just go with a Size for absolute safety (max upper bound on number rays traced can be assumed to be the number of points)
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>> real_pt_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> curr_intensity_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> next_intensity_;
    // pc::multi_threading::ThreadPrivate<Vector<Size>> freq_matching_;//for (mis)matching the frequency indices; wait: this assumes that we only go one positional step at a time, drastically increasing the difficulty for adding boundary conditions
    //NOTE: boundary conditions might result in
    //Again, no: this is merely for determining what freqs match, the actual previous point to use has not yet been defined!
    pc::multi_threading::ThreadPrivate<Matrix<Vector<Size>>> start_indices_;//For every point (on the ray) (maybe except the first point)), for every frequency, denotes from which index one needs to start
    // This is necessary due to boundary conditions possibly refererring to points way back on the ray; NOTE: vectors have size 2
    //TOOD: next thing not used anymore
    // pc::multi_threading::ThreadPrivate<Vector<unsigned char>> is_bdy_freq_;//for denoting which freqs (for the next point) should just be filled in with bdy conditions
    //either on a point-by-point basis (-> less memory), and is depending on the exact point either way (not used for anything but adding bdy conditions to the deque)
    //or we could compute it for all points at the same time
    //FIXME: for extreme overlap (a line quad completely covering another one), this definition is not well defined. TODO FIX THIS
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>> line_quad_left_overlap_;//for every line, denotes whether the left boundary overlaps with another quadrature
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>> line_quad_right_overlap_;//for denoting which freqs (for the next point) should just be filled in with bdy conditions
    //WAIT A MINUTE, as we only need the overlap in one direction (depending on the discretization direction),
    //NO, this is wrong, is there exist overlap in both directions; but this does not really matter that much, as this should only be used for the freq der boundary points
    //THUS TODO: replace with single thing
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>> line_quad_discdir_overlap_;//for denoting which lines should include the freq dir boundary conditions


    // If one precomputes S, Δτ, then one can set the boundary conditions by cheating (Δτ→very high, S→boundary intensity)
    // pc::multi_threading::ThreadPrivate<Vector<unsigned char>> bdy_freqs_;//for denoting which freqs (for the next point) should just be filled in with bdy conditions

    pc::multi_threading::ThreadPrivate<Matrix<Real>> intensities_;//for every point on the ray, for every frequency; stores the intensities for a single ray
    pc::multi_threading::ThreadPrivate<Matrix<Real>> bdy_intensities_;//for every point, for every split (two freqs per split); stores the intensities for a single ray
    pc::multi_threading::ThreadPrivate<Matrix<Real>> delta_tau_;//for every point, for every frequency, storing the optical depth increments; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>> S_;//for every point, for every frequency, storing source function; might be fiddled with to set boundary conditions
    //TODO: split S into S_curr and S_next, in order to freely change stuff for boundary conditions

    //For computing the second order accurate frequency derivative
    //Cant we just compute is during the main computation? I see no reason to store it now?
    //For the explicit part, this can certainly be done; however for the implicit part... (well it can also be computed on the fly, thus it might be useless)
    //Err no, it is necessary for allowing us to interpolate boundary values (using only 2 nonzero coeffs, third index may be arbitrary?, as third coef is arbitrary; DOCS STATE NO OOB ACCESS, so arbitrary in bounds value is fine)
    // Thus we might just shift the coef indices to set everything in a sensible range for the boundary interpolation stuff
    //These ones cannot be computed on the fly, as we kind of need to cheat including the boundary conditions efficiently
    // Weird enough, our boundary conditions can be written as special cases of the comoving RTE.
    // For initial boundary stuff, set coeffs to 0, Δτ→very large, S→bdy intensity
    // For interpolating previous values (linearly/cubibcly), set the coeffs respectively to the interpolation weights (Δτ should be 0)

    /// Explicit part (technically we should not need the coefs, however if we do not include this, we will access oob array values (even though the intention is to mutliply them by 0!))
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef1_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef2_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef3_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    ///
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index1_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index2_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index3_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions

    /// Implicit part (same remark)
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef1_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef2_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>> dIdnu_coef3_;//storing freq derivative coef; might be fiddled with to set boundary conditions
    ///
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index1_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index2_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>> dIdnu_index3_;//storing freq derivative corresp index; might be fiddled with to set boundary conditions
    // pc::multi_threading::ThreadPrivate<std::map




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
              Real&  chi ) const;

    template <ApproximationType approx>
    inline void compute_S_dtau_line_integrated (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext);

    inline Real compute_dtau_single_line(Model& model, Size curridx, Size nextidx, Size lineidx, Real curr_freq, Real next_freq, Real dz);

    template<ApproximationType approx>
    accel inline void compute_source_dtau (Model& model, Size currpoint, Size nextpoint, Size line, Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ, bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr, Real& Snext);


    accel inline void update_Lambda (
              Model &model,
        const Size   rr,
        const Size   f  );


    accel inline void solve_shortchar_order_0 (Model& model);
    accel inline void solve_shortchar_order_0 (
              Model& model,
        const Size   o,
        const Size   r);


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
    accel inline void set_eta_and_chi        (Model& model, const Size rr) const;
    accel inline void set_boundary_condition (Model& model) const;


    // Solvers for column densities
    ///////////////////////////////
    accel inline void set_column (Model& model) const;
    accel inline Real get_column (const Model& model, const Size o, const Size r) const;
};


#include "solver.tpp"
