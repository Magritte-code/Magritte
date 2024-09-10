#pragma once

#include "model/model.hpp"
#include "tools/types.hpp"

#include <map>
// #include <deque>//TODO: check if I actually use this thing anymore
#include <tuple>

///  Approximation used in the solver
/////////////////////////////////////
enum ApproximationType { None, OneLine, CloseLines };

struct Solver {

    // This seems to be the best value for the imager to avoid
    // dead pixels, while somewhat decently handling the inner
    // boundary
    // TODO: for the imager: let the ray stop at the inner bdy
    const Size MAX_CONSECUTIVE_BDY = 5; // for the new imager, we need some stopping
                                        // criterion to determine when the ray ends

    pc::multi_threading::ThreadPrivate<Vector<double>> dZ_; ///< distance increments along the ray
    pc::multi_threading::ThreadPrivate<Vector<Size>> nr_; ///< corresponding point number on the ray
    pc::multi_threading::ThreadPrivate<Vector<double>> shift_; ///< Doppler shift along the ray

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
    pc::multi_threading::ThreadPrivate<Vector<Real>> intensity_;

    // Comoving approach: TODO: clean up

    //TODO: make sure that the index ordering is consistent
    // For tracing the rays and keeping track of the closest ones
    Matrix<std::tuple<Size, Size>> corresponding_ray; // original point index, original direction index mapping to ray origin and the corresponding direction index
    Vector<std::map<std::tuple<Size, Size>, Size>> intensity_origin; // for every point, contains a map with as keys a set of origin points and corresponding directions (for the origins) and as values the direction index to be used for that point
    Vector<Vector<Size>> points_to_trace_ray_through; // raydir, rays to trace
    Matrix<Size>
        n_rays_through_point; // raydir, point ; might not be used that much in the future, as we
                              // might instead just use the closest ray for determining the result
    Matrix<Size> elements_in_rays_starting_from_origin; // raydir, point
    Vector<std::tuple<Size, Size>> rays_single_datapoint; // arbitrary idx to raydir, pointidx
    // hmm, we need to weight the rays somehow; for now we just only use the closest ray to
    // determine intensity
    // Matrix<Size> closest_ray;     // raydir, point (contains closest rayid?)
    Matrix<Size> min_ray_distsqr; // raydir, pointid
    // I will assume no more than 2^16-1 rays go through a specific point if choosing unsigned int
    // However, just go with a Size for absolute safety (max upper bound on number rays traced can
    // be assumed to be the number of points)
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>> real_pt_;

    /// Specific part for non-approximate comoving solver
    // pc::multi_threading::ThreadPrivate<Vector<Real>> curr_intensity_;
    // pc::multi_threading::ThreadPrivate<Vector<Real>> next_intensity_;
    // pc::multi_threading::ThreadPrivate<Vector<Size>> freq_matching_;//for (mis)matching the
    // frequency indices; wait: this assumes that we only go one positional step at a time,
    // drastically increasing the difficulty for adding boundary conditions
    // NOTE: boundary conditions might result in
    // Again, no: this is merely for determining what freqs match, the actual previous point to use
    // has not yet been defined!
    pc::multi_threading::ThreadPrivate<Matrix<Vector<Size>>>
        start_indices_; // For every point (on the ray) (maybe except the first point)), for every
                        // frequency, denotes from which index (pointonray, freq) one needs to start
                        // the computation
    // pc::multi_threading::ThreadPrivate<Matrix<std::array<Size, 2>>> start_indices_;//For every
    // point (on the ray) (maybe except the first point)), for every frequency, denotes from which
    // index (pointonray, freq) one needs to start the computation This is necessary due to boundary
    // conditions possibly refererring to points way back on the ray; NOTE: vectors have size 2
    // Done: use std::array<Size,2> instead of variable size Vector; errn never mind, need paracabs
    // Vector to properly map to gpu

    // WAIT A MINUTE, as we only need the overlap in one direction (depending on the discretization
    // direction), NO, this is wrong, is there exist overlap in both directions; but this does not
    // really matter that much, as this should only be used for the freq der boundary points THUS
    // TODO: replace with single thing To figure out which lines overlap with eachother (unsigned
    // char due to vector<bool> being a non-standard datatype)
    pc::multi_threading::ThreadPrivate<Vector<unsigned char>>
        line_quad_discdir_overlap_; // for denoting which lines should include the freq dir boundary
                                    // conditions

    // The overlap ranges
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        left_bound_; // Specifies the left bounds of the ranges in [Hz]
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        right_bound_; // Specifies the right bounds of the ranges in [Hz]
    pc::multi_threading::ThreadPrivate<Size> nb_ranges_; // contains the number of ranges
    pc::multi_threading::ThreadPrivate<Vector<Size>>
        left_bound_index_; // Specifies the corresponding freq index to the left bounds of the
                           // ranges
    pc::multi_threading::ThreadPrivate<Vector<Size>>
        right_bound_index_; // Specifies the corresponding freq index to the right bounds of the
                            // ranges

    // For computing the overlap ranges
    // TODO: use sorted lines to compute this instead: for two successive lines, check distance
    // compared to max line width; if already large enough, do not bother checking more lines;
    // otherwise check next line (this i due to the fact that a wide line might be 'hiding' behind a
    // peaked line)
    pc::multi_threading::ThreadPrivate<Vector<Size>>
        line_count_; // for every line, counts the number of quadrature points encountered
    pc::multi_threading::ThreadPrivate<Vector<Size>>
        quad_range_weight_; // for every quad, contains one if the quadrature may be counted for the
                            // range
    pc::multi_threading::ThreadPrivate<Size> tot_quad_range_weight_; // sum of quad_range_weight
    // TODO: for more general quadratures, change this
    // Err no, check outermost freqs of successive line in correct direction.
    // To make sure no 3 overlapping lines exist (with different line width), also check if the
    // lines lie far enough from eachother to be non-overlapping (use same criterion as usual?)
    //  pc::multi_threading::ThreadPrivate<Real> rootdiff_;//... in units of line width

    // If one precomputes S, Δτ, then one can set the boundary conditions by cheating (Δτ→very high,
    // S→boundary intensity) pc::multi_threading::ThreadPrivate<Vector<unsigned char>>
    // bdy_freqs_;//for denoting which freqs (for the next point) should just be filled in with bdy
    // conditions

    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        intensities_; // for every point on the ray, for every frequency; stores the intensities for
                      // a single ray
    // pc::multi_threading::ThreadPrivate<Matrix<Real>> bdy_intensities_;//for every point, for
    // every split (two freqs per split); stores the intensities for a single ray
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        delta_tau_; // for every point, for every frequency, storing the optical depth increments;
                    // might be fiddled with to set boundary conditions
    // pc::multi_threading::ThreadPrivate<Matrix<Real>> S_;//for every point, for every frequency,
    // storing sourcve function; might be fiddled with to set boundary conditions
    // TODO: split S into S_curr and S_next, in order to freely change stuff for boundary conditions
    // Seemingly strange to be duplicate, but this allows us to perfectly define boundary conditions
    // for each individual point
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        S_curr_; // for every point, for every frequency, storing source function; might be fiddled
                 // with to set boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        S_next_; // for every point, for every frequency, storing source function; might be fiddled
                 // with to set boundary conditions

    // For containing the boundary conditions more efficiently (for each line, one can then define
    // some ranges of boundary frequencies to find overlap faster)
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Real>>>
    //  left_bdy_deque_frequencies_;//contains the boundary frequencies
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Size>>>
    //  left_bdy_deque_rayposidx_;//contains the boundary frequencies ray position indices
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Size>>>
    //  left_bdy_deque_freq_idx_;//contains the boundary frequencies frequency indices
    //
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Real>>>
    //  right_bdy_deque_frequencies_;//contains the boundary frequencies
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Size>>>
    //  right_bdy_deque_rayposidx_;//contains the boundary frequencies ray position indices
    //  pc::multi_threading::ThreadPrivate<Vector<std::deque<Size>>>
    //  right_bdy_deque_freq_idx_;//contains the boundary frequencies frequency indices

    // For computing the second order accurate frequency derivative
    //  This needs to be done a priori, even though for default computation (no boundary
    //  conditions), it could be done on the fly. The main reason is that we can (mis)use these
    //  coefficients to construct boundary conditions into the same formulation as any regular
    //  (non-bdy) computation. Weirdly enough, our boundary conditions can be written as special
    //  cases of the comoving RTE. For initial boundary stuff, set coeffs to 0, Δτ→very large, S→bdy
    //  intensity For interpolating previous values (linearly/cubibcly), set the coeffs respectively
    //  to the interpolation weights (Δτ should be ≃0)

    /// Explicit part (technically we should not need the coefs, however if we do not include this,
    /// we will access oob array values (even though the intention is to mutliply them by 0!))
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef1_curr_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef2_curr_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef3_curr_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    /// Note: these indices contain only the FREQUENCY index, NOT THE POINT INDEX, (as that should
    /// be derived from start_indices_)
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index1_curr_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index2_curr_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index3_curr_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions

    /// Implicit part (same remark)
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef1_next_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef2_next_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    pc::multi_threading::ThreadPrivate<Matrix<Real>>
        dIdnu_coef3_next_; // storing freq derivative coef; might be fiddled with to set boundary
                           // conditions
    /// Note: these indices contain only the FREQUENCY index, NOT THE POINT INDEX, (as that is
    /// pointless in this case (same nextpoint as one need to fill into this matrix)
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index1_next_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index2_next_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions
    pc::multi_threading::ThreadPrivate<Matrix<Size>>
        dIdnu_index3_next_; // storing freq derivative corresp index; might be fiddled with to set
                            // boundary conditions
    // EXTRA REQUIREMENT: dIdnu_index1_next_()==start_indices_()[1]; this for easily treating the
    // implicit part
    //  && dIdnu_index2/3_next_()!=start_indices_()[1] IF the corresponding coefficient is nonzero
    // Practically, this means that I can just simply subtract

    /// Specific part for comoving approximate solver (CoMoving Approximate)
    // As this only deals with a single position increment at a time, less data needs to be stored
    // compared to the non-approximate version unless stated otherwise, all are indexed using the
    // sorted frequency index; however, which sorted index (curr or next point) depends on what is
    // stored TODO EXPLAIN
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_computed_intensities_; // stores the computed intensities//TODO: INIT WITH CMB
    // pc::multi_threading::ThreadPrivate<Vector<Real>> cma_curr_frequencies_;//stores the
    // frequencies of the current point
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_start_intensities_; // stores the intensity at the start of increment
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_start_frequencies_; // stores the start frequency (mainly convenient for programming)
    pc::multi_threading::ThreadPrivate<Vector<Size>>
        cma_start_frequency_index_; // stores the starting frequency index (used for looking up
                                    // stuff)
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_end_frequencies_; // stores the end frequency (mainly convenient for programming)
    // pc::multi_threading::ThreadPrivate<Vector<Real>> cma_end_frequency_index_;//stores the end
    // frequency index (used for looking up stuff)
    // not necessary to store, as this maps the end frequency index onto itself.
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_chi_curr_; // stores the starting opacity (thus does not need to be computed twice)
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_chi_next_; // stores the ending opacity. Needs to overwrite cma_chi_curr when the data
                       // is fully mapped
    pc::multi_threading::ThreadPrivate<Vector<char>>
        cma_compute_curr_opacity_; // stores whether to compute the opacity
    pc::multi_threading::ThreadPrivate<Vector<char>>
        cma_compute_next_opacity_; // stores whether to compute the opacity. Needs to overwrite
                                   // cma_compute_curr_opacity_ when data is fully mapped
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_S_curr_; // stores the computed current source function for the increment
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_S_next_; // stores the computed next source function for the increment
    // TODO: remove these coef storing things
    //  pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef1_curr_;//stores derivative
    //  coefficents for the intensity derivative with respect to the frequency
    //  pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef2_curr_;//1 stores the coef
    //  of the current freq, 2 of the nearby freq and 3 of the outermost freq
    //  pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef3_curr_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef1_next_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef2_next_;
    pc::multi_threading::ThreadPrivate<Vector<Real>> cma_dIdnu_coef3_next_;
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_dIdnu_expl_; // stores the computed explicit frequency derivative
    pc::multi_threading::ThreadPrivate<Vector<Real>>
        cma_delta_tau_; // stores the optical depth increments belonging to each position increment

    Vector<Real> eta;
    Vector<Real> chi;

    Size nblocks  = 512;
    Size nthreads = 512;

    Size length;
    Size centre;
    Size width;

    Size n_off_diag;

    template <Frame frame, bool use_adaptive_directions> void setup(Model& model);

    // template <Frame frame>
    void setup_new_imager(Model& model, Image& image, const Vector3D& ray_dir);
    void setup(const Size l, const Size w, const Size n_o_d);

    template <bool use_adaptive_directions> void setup_comoving(Model& model);
    void setup_comoving_new_imager(Model& model, Image& image, const Vector3D& ray_dir);
    void setup_comoving(Model& model, const Size length, const Size width);

    accel inline Real get_dshift_max(const Model& model, const Size o);

    template <Frame frame, bool use_adaptive_directions> inline void get_ray_lengths(Model& model);

    template <Frame frame, bool use_adaptive_directions>
    inline Size get_ray_lengths_max(Model& model);

    // template <Frame frame>
    inline Size get_ray_lengths_max_new_imager(Model& model, Image& image, const Vector3D& ray_dir);

    accel inline Size get_ray_length_new_imager(const Geometry& geometry, const Vector3D& origin,
        const Size start_bdy, const Vector3D& raydir);

    template <Frame frame, bool use_adaptive_directions>
    accel inline Size trace_ray(const Geometry& geometry, const Size o, const Size r,
        const double dshift_max, const int increment, Size id1, Size id2);

    // With extra functionality to figure out when to stop our computations on the ray
    template <Frame frame, bool use_adaptive_directions>
    accel inline Size trace_ray_comoving(const Geometry& geometry, const Size o, const Size r,
        const Size rr, const Size rayidx, const double dshift_max, const int increment, Size id1,
        Size id2, Size& outermost_interesting_point_rayidx);

    accel inline Size trace_ray_imaging_get_start(const Geometry& geometry, const Vector3D& origin,
        const Size start_bdy, const Vector3D& raydir, Real& Z);

    template <Frame frame>
    accel inline Size trace_ray_imaging(const Geometry& geometry, const Vector3D& origin,
        const Size start_bdy, const Vector3D& raydir, const double dshift_max, const int increment,
        Real& Z, Size id1, Size id2);

    accel inline void set_data(const Size crt, const Size nxt, const double shift_crt,
        const double shift_nxt, const double dZ_loc, const double dshift_max, const int increment,
        Size& id1, Size& id2);

    accel inline Real gaussian(const Real width, const Real diff) const;
    accel inline Real planck(const Real temp, const Real freq) const;

    accel inline Real boundary_intensity(const Model& model, const Size p, const Real freq) const;

    template <ApproximationType approx>
    accel inline void get_eta_and_chi(const Model& model, const Size p, const Size l,
        const Real freq, Real& eta, Real& chi) const;

    template <ApproximationType approx>
    inline void compute_S_dtau_line_integrated(Model& model, Size currpoint, Size nextpoint,
        Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext);

    inline Real compute_dtau_single_line(Model& model, Size curridx, Size nextidx, Size lineidx,
        Real curr_freq, Real next_freq, Real dz);

    template <ApproximationType approx>
    accel inline void compute_source_dtau(Model& model, Size currpoint, Size nextpoint, Size line,
        Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ,
        bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr,
        Real& Snext);

    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void update_Lambda(Model& model, const Size rr, const Size f);
    // accel inline void solve_shortchar_order_0_ray_forward (
    //           Model& model,
    //           const Size   o,
    //           const Size   r);
    // accel inline void solve_shortchar_order_0_ray_backward
    // (
    //           Model& model,
    //           const Size   o,
    //           const Size   r);

    // Comoving solvers stuff
    /////////////////////////
    // general ray tracing differences
    //  inline void setup_comoving (Model& model, const Size l, const Size w);
    template <bool use_adaptive_directions>
    inline void get_static_rays_to_trace(Model& model);
    template <bool use_adaptive_directions>
    accel inline void trace_ray_points(const Geometry& geometry, const Size o, const Size rdir,
        const Size rsav, const Size rayidx);
    // Complicated solver stuff

    inline void match_frequency_indices(Model& model, const Size nextpoint, const Size currpoint,
        const Real next_shift, const Real curr_shift, Size nextpointonrayindex,
        Size currpointonrayindex, bool is_upward_disc);
    inline void get_overlapping_lines(Model& model, const Size nextpoint, bool is_upward_disc);
    inline void get_line_ranges(
        Model& model, const Size curr_point, bool is_upward_disc, Real curr_shift);
    inline void set_implicit_boundary_frequencies(Model& model, const Size nextpoint,
        const Size nextpointonrayindex, Real shift_next,
        std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index,
        bool is_upward_disc);
    inline void match_overlapping_boundary_conditions(Model& model, const Size currpoint,
        const Size curr_point_on_ray_index, const Real curr_shift,
        std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index);
    inline void set_initial_boundary_conditions(Model& model, const Size initial_bdy,
        const Real curr_shift,
        std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index);
    template <ApproximationType approx>
    inline void comoving_ray_bdy_setup_forward(Model& model, Size first_interesting_rayposidx);
    template <ApproximationType approx>
    inline void comoving_ray_bdy_setup_backward(Model& model, Size last_interesting_rayposidx);
    template <ApproximationType approx, bool use_adaptive_directions>
    inline void solve_comoving_order_2_sparse(Model& model);
    template <bool use_adaptive_directions>
    inline void solve_comoving_single_step(Model& model, const Size rayposidx, const Size o,
        const Size rr, const bool is_upward_disc, const bool forward_ray);
    template <ApproximationType approx, bool use_adaptive_directions>
    inline void solve_comoving_order_2_sparse(Model& model,
        const Size o,      // ray origin point
        const Size r,      // ray direction index
        const Size rayidx, // ray trace index
        const double dshift_max);
    accel inline Size trace_ray_indicate_point(const Geometry& geometry, const Size o, const Size r,
        const double dshift_max, const int increment, Size id1, Size id2);
    accel inline void set_data_indicate_point(const Size crt, const Size nxt,
        const double shift_crt, const double shift_nxt, const double dZ_loc,
        const double dshift_max, const int increment, Size& id1, Size& id2);

    // Comoving approx solver stuff
    template <bool use_adaptive_directions>
    accel inline void setup_comoving_local_approx(Model& model);
    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_comoving_local_approx_order_2_sparse(Model& model);
    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_comoving_local_approx_order_2_sparse(Model& model,
        const Size o,      // ray origin point
        const Size r,      // ray direction index
        const Size rayidx, // ray trace index
        const double dshift_max);
    template <ApproximationType approx>
    accel inline void comoving_local_approx_map_data(Model& model, const Size curr_point,
        const Size next_point, const Real curr_shift, const Real next_shift,
        const bool is_upward_disc, const Real dZ, const Size bdy_point);
    template <ApproximationType approx>
    accel inline void comoving_approx_map_single_data(Model& model, const Size curr_point,
        const Size next_point, const Real dZ, const bool is_in_bounds, const Size curr_freq_count,
        const Size next_freq_idx, const Size curr_freq_idx, const Real next_freq,
        const Real curr_freq, const Real next_shift, const Real curr_shift,
        const Size curr_line_idx, const bool is_upward_disc, const Size bdy_point);
    template <bool use_adaptive_directions>
    accel inline void solve_comoving_local_approx_single_step(Model& model, const Size next_point,
        const Size o, const Size rr, const bool is_upward_disc);

    // Point pruning solvers stuff
    //////////////////////////////

    accel inline bool check_close_line(const Real currfreq, const Real nextfreq,
        const Size currpoint, const Size nextpoint, const Model& model);
    // accel inline bool check_close_line (const Real prevfreq, const Real currfreq, const Real
    // nextfreq, const Size prevpoint, const Size currpoint, const Size nextpoint, const Model&
    // model);
    template <Frame frame, bool use_adaptive_directions>
    accel inline Size trace_ray_pruned(const Model& model, const Size o, const Size r,
        const double dshift_max, const int increment, Size id1, Size id2, const Real freq);

    template <ApproximationType approx, bool use_adaptive_directions>
    inline void solve_feautrier_order_2_sparse_pruned_rays(Model& model);

    // Solvers for images
    /////////////////////
    // algorithms for tracing the rays and extracting the
    // information of the solver
    template <ApproximationType approx>
    accel inline void image_feautrier_order_2(Model& model, const Size rr);
    template <ApproximationType approx>
    inline void image_feautrier_order_2_new_imager(
        Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix);
    template <ApproximationType approx>
    inline void image_shortchar_order_0_new_imager(
        Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix);
    // actual solver
    template <ApproximationType approx>
    accel inline void image_feautrier_order_2(Model& model, const Size o, const Size f);
    template <ApproximationType approx>
    accel inline Real image_shortchar_order_0(Model& model, const Size o, const Size f);

    template <ApproximationType approx>
    accel inline void image_feautrier_order_2_for_point(Model& model, const Size rr, const Size p);
    template <ApproximationType approx>
    accel inline void image_feautrier_order_2_for_point_loc(
        Model& model, const Size o, const Size f);

    template <ApproximationType approx>
    accel inline void image_optical_depth(Model& model, const Size rr);
    template <ApproximationType approx>
    inline void image_optical_depth_new_imager(
        Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix);
    // actual solver
    template <ApproximationType approx>
    accel inline void image_optical_depth(Model& model, const Size o, const Size f);

    // Comoving solver
    //////////////////
    template <ApproximationType approx>
    inline void image_comoving_new_imager(Model& model, const Vector3D& ray_dir, const Size nxpix,
        const Size nypix, const Vector<Real>& image_freqs);

    // actual solver; assumes the full ray computation has been set up correctly
    inline void solve_comoving_image_single_step(
        Model& model, const Size rayposidx, const bool is_upward_disc, const bool forward_ray);

    // interpolation of computed comoving intensities
    inline void interpolate_computed_comoving_intensities(Model& model, Image& image, Size pixidx,
        const Size currpoint, const Size curr_point_on_ray_index, const Real curr_shift,
        std::multimap<Real, Size>& multimap_image_freq_to_index);

    // Solvers only computing u
    ///////////////////////////
    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_feautrier_order_2(Model& model);

    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_feautrier_order_2_sparse(Model& model);

    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_feautrier_order_2_anis(Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2(Model& model, const Size o, const Size f);

    // // Solvers for both u and v
    // ///////////////////////////

    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_shortchar_order_0(Model& model);
    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_shortchar_order_0(Model& model, const Size o, const Size r);
    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_shortchar_order_0_sparse(Model& model, const Size o, const Size r);

    template <ApproximationType approx, bool use_adaptive_directions>
    accel inline void solve_feautrier_order_2_uv(Model& model);

    template <ApproximationType approx>
    accel inline void solve_feautrier_order_2_uv(Model& model, const Size o, const Size f);

    // Getters for emissivities, opacities, and boundary
    // conditions
    ///////////////////////////////////////////////////////////////
    accel inline void set_eta_and_chi(Model& model, const Size rr) const;
    accel inline void set_boundary_condition(Model& model) const;

    // Solvers for column densities
    ///////////////////////////////
    template <bool use_adaptive_directions> accel inline void set_column(Model& model) const;
    template <bool use_adaptive_directions>
    accel inline Real get_column(const Model& model, const Size o, const Size r) const;
};

#include "solver.tpp"
