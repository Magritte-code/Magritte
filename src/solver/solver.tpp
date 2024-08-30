#include <set>
#include <unordered_map>

///  Prepares datastructures for NLTE calculations, old imager using regular solvers (shortchar,
///  feautrier)
///   @param[in] model: model to apply the solver to
template <Frame frame, bool use_adaptive_directions> inline void Solver ::setup(Model& model) {
    const Size length = 2 * get_ray_lengths_max<frame, use_adaptive_directions>(model) + 1;
    const Size width  = model.parameters->nfreqs();
    const Size n_o_d  = model.parameters->n_off_diag;

    model.set_dshift_max(); // err, probably belongs somewhere
                            // else, but we need to compute
                            // the max shift for each point

    setup(length, width, n_o_d);
}

///  Prepares datastructures for the new imager using the regular solvers (shortchar, feautrier)
///   @param[in] model: model to image
///   @param[in] image: image object to store data into
///   @param[in] ray_dir: ray direction for image
inline void Solver ::setup_new_imager(Model& model, Image& image, const Vector3D& ray_dir) {
    const Size length = 2 * get_ray_lengths_max_new_imager(model, image, ray_dir) + 1;
    const Size width  = model.parameters->nfreqs();
    const Size n_o_d  = model.parameters->n_off_diag;

    model.set_dshift_max(); // err, probably belongs somewhere
                            // else, but we need to compute
                            // the max shift for each point

    setup(length, width, n_o_d);
}

///  Internal helper method for preparing datastructures for regular solvers (shortchar, feautrier)
///   @param[in] l: "length": maximum ray length*2+1; this ensures that all indices pertaining to
///   position on a ray are in bounds
///   @param[in] w: "width": number of frequencies; ensures that all indices pertaining to
///   frequencies are in bounds
///   @param[in] n_o_d: number of off-diagonal elements; no longer used
inline void Solver ::setup(const Size l, const Size w, const Size n_o_d) {
    length     = l;
    centre     = l / 2;
    width      = w;
    n_off_diag = n_o_d;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++) {
        // general data structures, NOTE: not all datastructures are used
        dZ_(i).resize(length);    // distance increment
        nr_(i).resize(length);    // indices of point
        shift_(i).resize(length); // doppler shift

        eta_c_(i).resize(width); // emissivity at current point
        eta_n_(i).resize(width); // emissivity at next point

        chi_c_(i).resize(width); // opacity at current point
        chi_n_(i).resize(width); // opacity at next point

        source_c_(i).resize(width); // source function at current point
        source_n_(i).resize(width); // source function at next point

        inverse_chi_(i).resize(length); // 1/opacity

        tau_(i).resize(width); // total optical depth

        // Helper variables for solving feautrier
        Su_(i).resize(length);
        Sv_(i).resize(length);

        A_(i).resize(length);
        C_(i).resize(length);
        inverse_A_(i).resize(length);
        inverse_C_(i).resize(length);

        FF_(i).resize(length);
        FI_(i).resize(length);
        GG_(i).resize(length);
        GI_(i).resize(length);
        GP_(i).resize(length);

        // Accelerated lambda iteration helper variables
        L_diag_(i).resize(length);

        L_upper_(i).resize(n_off_diag, length);
        L_lower_(i).resize(n_off_diag, length);
    }
}

///  Prepares datastructures for NLTE calculations using the comoving solver
///   @param[in] model: model to apply the solver to
template <bool use_adaptive_directions>
inline void Solver::setup_comoving(Model& model) {
    length = 2 * get_ray_lengths_max<Rest, use_adaptive_directions>(model) + 1;
    width  = model.parameters->nfreqs();
    // const Size n_o_d = model.parameters->n_off_diag;

    model.set_dshift_max(); // err, probably belongs somewhere else, but we need to compute the max
    // shift for each point

    setup_comoving(model, length, width);

    // Finally also trace the rays in advance to prune the unnecessary ones.
    get_static_rays_to_trace<use_adaptive_directions>(model);
}

///  Prepares datastructures for the new imager using the comoving solver
///   @param[in] model: model to image
///   @param[in] image: image object to store data into
///   @param[in] ray_dir: ray direction for image
inline void Solver::setup_comoving_new_imager(Model& model, Image& image, const Vector3D& ray_dir) {
    length = 2 * get_ray_lengths_max_new_imager(model, image, ray_dir) + 1;
    width  = model.parameters->nfreqs();
    // const Size n_o_d  = model.parameters->n_off_diag; offdiagonal accelerated
    // lambda iteration is no longer used and is too complicated for the comoving solver

    model.set_dshift_max(); // err, probably belongs somewhere
                            // else, but we need to compute
                            // the max shift for each point

    setup_comoving(model, length, width);
}

///  Internal helper method for preparing datastructures for the comoving solver
///   @param[in] length: maximum ray length*2+1; this ensures that all indices pertaining to
///   position on a ray are in bounds
///   @param[in] width: number of frequencies; ensures that all indices pertaining to frequencies
///   are in bounds
inline void Solver ::setup_comoving(Model& model, const Size length, const Size width) {
    centre = length / 2; // length, width should be already setup in the function calling this one

    model.set_dshift_max(); // err, probably belongs somewhere else, but we need to compute the max
                            // shift for each point

    points_to_trace_ray_through.resize(model.parameters->hnrays());
    for (Size i = 0; i < model.parameters->hnrays(); i++) {
        points_to_trace_ray_through[i].resize(model.parameters->npoints());
    }
    // TODO: implement new datastructures for interpolating the computed intensities (healpix based)
    // Also implement the functions for getting the corresponding directions for each individual point (if it can be mapped) 
    corresponding_ray.resize(model.parameters->npoints(), model.parameters->hnrays());
    intensity_origin.resize(model.parameters->npoints());



    // For determining which ray lies closest to each point
    n_rays_through_point.resize(model.parameters->hnrays(), model.parameters->npoints());
    min_ray_distsqr.resize(model.parameters->hnrays(), model.parameters->npoints());
    // closest_ray.resize(model.parameters->hnrays(), model.parameters->npoints());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++) {
        // general ray tracing variables
        dZ_(i).resize(length);    // distance increment
        nr_(i).resize(length);    // point index
        shift_(i).resize(length); // doppler shift

        start_indices_(i).resize(length, width);
        for (Size j = 0; j < length; j++) {
            for (Size k = 0; k < width; k++) {
                start_indices_(i)(j, k).resize(
                    2); // for every thread i, for every ray position j, frequency index k, contains
                        // a tuple with the starting ray position, frequency index
            }
        }

        line_quad_discdir_overlap_(i).resize(
            model.parameters->nlines()); // whether the line quadratures overlap
        left_bound_(i).resize(
            model.parameters->nlines()); // left frequency bounds (for each continguous line region)
                                         // of the current frequency discretization
        right_bound_(i).resize(
            model.parameters->nlines()); // right frequency bounds (for each continguous line
                                         // region) of the current frequency discretization
        left_bound_index_(i).resize(
            model.parameters->nlines()); // corresponding frequency index for left bound
        right_bound_index_(i).resize(
            model.parameters->nlines()); // corresponding frequency index for right bound

        line_count_(i).resize(model.parameters->nlines()); // helper variable for counting number of
                                                           // quadratures encountered for each line
        quad_range_weight_(i).resize(
            model.parameters->nquads()); // helper variable for counting number of
                                         // quadratures encountered for each line; contains weight 1
                                         // for non-bdy frequencies, 0 for boundary frequencies

        intensities_(i).resize(
            length, width); // stores computed intensities at each ray position, frequency index
        delta_tau_(i).resize(length, width); // stores computed optical depth increments at each ray
                                             // position, frequency index
        S_curr_(i).resize(length, width);    // stores the current source function
        S_next_(i).resize(length, width);    // stores the next source function

        // Stores coefficients for computing the intensity derivatives with respect to the frequency
        dIdnu_coef1_curr_(i).resize(length, width);
        dIdnu_coef2_curr_(i).resize(length, width);
        dIdnu_coef3_curr_(i).resize(length, width);

        dIdnu_index1_curr_(i).resize(length, width);
        dIdnu_index2_curr_(i).resize(length, width);
        dIdnu_index3_curr_(i).resize(length, width);

        dIdnu_coef1_next_(i).resize(length, width);
        dIdnu_coef2_next_(i).resize(length, width);
        dIdnu_coef3_next_(i).resize(length, width);

        dIdnu_index1_next_(i).resize(length, width);
        dIdnu_index2_next_(i).resize(length, width);
        dIdnu_index3_next_(i).resize(length, width);
    }
}

///  Prepares datastructures for NLTE calculations for the approximate comoving solver
///   @param[in] model: model to apply the solver to
template <bool use_adaptive_directions>
inline void Solver ::setup_comoving_local_approx(Model& model) {
    length = 2 * get_ray_lengths_max<Rest, use_adaptive_directions>(model) + 1;
    centre = length / 2;
    width  = model.parameters->nfreqs();

    model.set_dshift_max(); // err, probably belongs somewhere else, but we need to compute the max
                            // shift for each point

    std::cout << "setup length: " << length << std::endl;
    std::cout << "setup width: " << width << std::endl;

    // Setup for allowing a ray to compute multiple values
    points_to_trace_ray_through.resize(model.parameters->hnrays());
    for (Size i = 0; i < model.parameters->hnrays(); i++) {
        points_to_trace_ray_through[i].resize(model.parameters->npoints());
    }
    corresponding_ray.resize(model.parameters->npoints(), model.parameters->hnrays());
    intensity_origin.resize(model.parameters->npoints());

    // for determining which ray lies closest to each point
    n_rays_through_point.resize(model.parameters->hnrays(), model.parameters->npoints());
    min_ray_distsqr.resize(model.parameters->hnrays(), model.parameters->npoints());
    // closest_ray.resize(model.parameters->hnrays(), model.parameters->npoints());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++) {
        // general ray tracing variables
        dZ_(i).resize(length);    // distance increment
        nr_(i).resize(length);    // point index
        shift_(i).resize(length); // doppler shift

        // For simplicity, we define other variables for the CoMoving Approx solver, as we store
        // less (in the full solver, we use matrices, here vectors suffice)
        cma_computed_intensities_(i).resize(width); // contains the currently computed intensities
        cma_start_intensities_(i).resize(width);    // contains the previously computed intensities
        cma_start_frequencies_(i).resize(
            width); // contains the rest frame frequencies at the current point
        cma_start_frequency_index_(i).resize(
            width); // contains the frequency indices for the comoving solver at the current point
        cma_end_frequencies_(i).resize(
            width);                     // contains the rest frame frequencies at the next point
        cma_chi_curr_(i).resize(width); // contains the opacities at the current point
        cma_chi_next_(i).resize(width); // contains the opacities at the next point
        cma_S_curr_(i).resize(width);   // contains the source functions at the current point
        cma_S_next_(i).resize(width);   // contains the source functions at the next point

        // contains bools: whether to recompute the opacities; necessary for compute_source_dtau
        cma_compute_curr_opacity_(i).resize(width);
        cma_compute_next_opacity_(i).resize(width);

        // helper variables for computing the intensity derivatives with respect to the frequency
        cma_dIdnu_coef1_next_(i).resize(width);
        cma_dIdnu_coef2_next_(i).resize(width);
        cma_dIdnu_coef3_next_(i).resize(width);
        cma_dIdnu_expl_(i).resize(width);
        cma_delta_tau_(i).resize(width);
    }

    // Finally also trace the rays in advance to prune the unnecessary ones.
    get_static_rays_to_trace<use_adaptive_directions>(model);
}

///  Getter for the maximum allowed shift value determined by the smallest line
///    @param[in] o : number of point under consideration
///    @retrun maximum allowed shift value determined by the smallest line
///////////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::get_dshift_max(const Model& model, const Size o) {
    Real dshift_max = std::numeric_limits<Real>::max();

    for (const LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        const Real inverse_mass   = lspec.linedata.inverse_mass;
        const Real new_dshift_max = model.parameters->max_width_fraction
                                  * model.thermodynamics.profile_width(inverse_mass, o);

        if (dshift_max > new_dshift_max) {
            dshift_max = new_dshift_max;
        }
    }

    return dshift_max;
}

///  Computes the maximum ray lengths when raytracing in a model
///    @param[in] model : the model
///////////////////////////////////////////////////////////////////////////////
template <Frame frame, bool use_adaptive_directions>
inline void Solver ::get_ray_lengths(Model& model) {
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar = model.geometry.rays.get_antipod_index(rr);

            const Real dshift_max = get_dshift_max(model, o);

            model.geometry.lengths(rr, o) =
                model.geometry.get_ray_length<frame, use_adaptive_directions>(o, rr, dshift_max)
                + model.geometry.get_ray_length<frame, use_adaptive_directions>(o, ar, dshift_max);
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec();
}

///  Getter for the maximum ray length needed for ray-tracing in a model
///   @param[in] model: the model
///   @returns the maximum ray length
template <Frame frame, bool use_adaptive_directions>
inline Size Solver ::get_ray_lengths_max(Model& model) {
    get_ray_lengths<frame, use_adaptive_directions>(model);

    Geometry& geo = model.geometry;

    geo.lengths_max = *std::max_element(geo.lengths.vec.begin(), geo.lengths.vec.end());

    return geo.lengths_max;
}

///  Traces all points on a given ray and checks whether the ray lies closest (among currently
///  traced rays) to the given ray. Should be called in both directions, but as tracing antipodal
///  rays should be symmetric, we only need to keep track of one direction. Used in setup for doing
///  NLTE computations using the (approximate) comoving solver
template <bool use_adaptive_directions>
accel inline void Solver ::trace_ray_points(const Geometry& geometry,
    const Size o,      // origin point of ray
    const Size rdir,   // ray direction to trace ∈ [0, nrays-1]
    const Size rsav,   // ray index save direction ∈ [0,hnrays-1]
    const Size rayidx) // for indexing the rays (counts from 0 to total number of rays traced
                       // through the domain
{
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = geometry.get_next<use_adaptive_directions>(o, rdir, o, Z, dZ);

    if (geometry.valid_point(nxt)) {
        // get distance and check if closest ray
        Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);
        // Compute the direction index at nxt (if valid).
        std::tuple<bool, Size> valid_rcur = geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rsav, nxt);
        bool valid = std::get<0>(valid_rcur);
        Size rcur  = std::get<1>(valid_rcur);
        if (valid)
        {
            // If it is the first time we encounter this point, or this is the closest ray: assign this
            // ray to compute the stuff (J, lambda) of the point
            if (n_rays_through_point(rcur, nxt) == 0 || dist2 < min_ray_distsqr(rcur, nxt)) {
                min_ray_distsqr(rcur, nxt) = dist2;
                corresponding_ray(nxt, rcur) = {o, rsav};
                // closest_ray(rcur, nxt)     = rayidx;
            }

            n_rays_through_point(rcur, nxt)++;
        }

        Size crt = o;

        while (geometry.not_on_boundary(nxt)) {
            crt = nxt;
            nxt = geometry.get_next<use_adaptive_directions>(o, rdir, nxt, Z, dZ);

            // get distance and check if closest ray
            Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);

            // Compute the direction index at nxt (if valid).
            std::tuple<bool, Size> valid_rcur = geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rsav, nxt);
            bool valid = std::get<0>(valid_rcur);
            Size rcur  = std::get<1>(valid_rcur);
            if (valid &&(n_rays_through_point(rsav, nxt) == 0 || dist2 < min_ray_distsqr(rsav, nxt))) {
                min_ray_distsqr(rsav, nxt) = dist2;
                corresponding_ray(nxt, rcur) = {o, rsav};
                // closest_ray(rsav, nxt)     = rayidx;
            }

            n_rays_through_point(rsav, nxt)++;
        }
    }
}

///  For all directions, determines a ray covering of the points for the (approximate) comoving
///  solver
///   @param[in] model: the model for which we want to use the comoving solver
template <bool use_adaptive_directions>
inline void Solver ::get_static_rays_to_trace(Model& model) {
    Vector<Size> n_points_to_trace_ray_through(model.parameters->hnrays());

    accelerated_for(rr, model.parameters->hnrays(), {
        Size n_rays_to_trace = 0;
        const Size ar        = model.geometry.rays.get_antipod_index(rr);
        // To make sure that we order the rays starting from the largest elements
        for (Size pointidx = 0; pointidx < model.parameters->npoints(); pointidx++) {
            const Size o = model.geometry.sorted_position_indices[pointidx];

            // DEBUG: commenting this out results in tracing through all points
            if (n_rays_through_point(rr, o) > 0) {
                continue;
            }

            // seemingly no ray has been traced through this point, so we must trace a ray through
            // it
            points_to_trace_ray_through[rr][n_rays_to_trace] = o;
            // tracing rays is symmetric, so only keep for the first half of the ray directions

            // trace ray through point
            n_rays_through_point(rr, o)++;
            // antipod has exactly same number of rays through the point, so do not save

            // For generality, I assign a ray index to each ray of a given direction rsav
            // Also the ray direction is identified with the forward dir
            // closest_ray(rr, o)     = n_rays_to_trace;
            corresponding_ray(o, rr) = {o, rr};

            min_ray_distsqr(rr, o) = 0.0;

            // now trace ray through rest of model
            trace_ray_points<use_adaptive_directions>(model.geometry, o, rr, rr, n_rays_to_trace);
            trace_ray_points<use_adaptive_directions>(model.geometry, o, ar, rr, n_rays_to_trace);

            n_rays_to_trace++;
        }
        n_points_to_trace_ray_through[rr]=n_rays_to_trace;
        points_to_trace_ray_through[rr].resize(
            n_rays_to_trace); // and now the correct size, instead of parameters->npoints()
    })

    // Convert the datastructure to something more usable
    for (Size idx = 0; idx<model.parameters->npoints(); idx++) {
        for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
            intensity_origin[idx][corresponding_ray(idx, rr)] = rr;
        }
    }

    // debug print stuff
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
        std::cout << "rr: " << rr
                  << " size points_to_trace_ray_through: " << points_to_trace_ray_through[rr].size()
                  << std::endl;
        for (Size idx = 0; idx < points_to_trace_ray_through[rr].size(); idx++) {
            // std::cout<<"point: "<<points_to_trace_ray_through[rr][idx]<<std::endl;
        }
        // std::cout<<"number of rays per point"<<std::endl;
        for (Size p = 0; p < model.parameters->npoints(); p++) {
            // std::cout<<"closest ray: "<<closest_ray(rr,p)<<std::endl;
            // std::cout<<"point: "<<p<<"#: "<<n_rays_through_point(rr, p)<<std::endl;
        }
    }
}

///  Comoving solver setup helper function: Matches the sorted frequency indices in a single
///  interval such that the frequency difference is minimal
///   @param[in] model: the model for which we want to apply the comoving solver
///   @param[in] nextpoint: the point index of the next point
///   @param[in] currpoint: the point index of the current point
///   @param[in] next_shift: doppler shift at the next point
///   @param[in] curr_shift: doppler shift at the current point
///   @param[in] nextpointonrayindex: the ray position index of the next point
///   @param[in] currpointonrayindex: the ray position index of the current point
///   @param[in] is_upward_disc: whether or not to use the rightmost frequency indices as boundary
///   conditions
inline void Solver ::match_frequency_indices(Model& model, const Size nextpoint,
    const Size currpoint, const Real next_shift, const Real curr_shift, Size nextpointonrayindex,
    Size currpointonrayindex, bool is_upward_disc) {
    // if upward discretization, we start from the uppermost part
    if (is_upward_disc) {
        // Starting from the highest frequency
        //+1 to all indices due to using unsigned ints when looping down (overflow otherwise)
        Size curr_freq_idx      = model.parameters->nfreqs() - 1 + 1;
        Size next_freq_idx      = model.parameters->nfreqs() - 1 + 1;
        const Size min_freq_idx = 0 + 1;
        // assumes at least a single frequency will not be a boundary condition (is reasonable if we
        // limit the doppler shift) NOTE: oob freqs will be dealt with later on
        while (curr_freq_idx >= min_freq_idx && next_freq_idx >= min_freq_idx) {
            // check if the frequency at the previous point index is higher than this frequency (in
            // static frame)
            if (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx - 1) * curr_shift
                > model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx - 1)
                      * next_shift) {
                curr_freq_idx--;
            } else // freq matching index must be determined if it is just higher than the freq at
                   // the previous point
            {
                start_indices_()(nextpointonrayindex, next_freq_idx - 1)[0] =
                    currpointonrayindex; // stores point on ray index
                start_indices_()(nextpointonrayindex, next_freq_idx - 1)[1] =
                    curr_freq_idx - 1; // stores freq
                // for the sake of boundary conditions, these start_indices_ might be overwritten
                // later
                next_freq_idx--;
            }
        }
        // match remaining points as best as possible (with the outermost point)
        while (next_freq_idx >= min_freq_idx) {
            start_indices_()(nextpointonrayindex, next_freq_idx - 1)[0] =
                currpointonrayindex; // stores point
            start_indices_()(nextpointonrayindex, next_freq_idx - 1)[1] =
                min_freq_idx - 1; // stores freq
            next_freq_idx--;
        }
    } else {
        // Starting from the lowest frequency
        Size curr_freq_idx      = 0;
        Size next_freq_idx      = 0;
        const Size max_freq_idx = model.parameters->nfreqs() - 1;
        // assumes at least a single frequency will not be a boundary condition (is reasonable if we
        // limit the doppler shift)
        while (curr_freq_idx <= max_freq_idx && next_freq_idx <= max_freq_idx) {
            // check if a previous frequency index exists which is higher than this (in static
            // frame)
            if (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx) * curr_shift
                < model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx) * next_shift) {
                curr_freq_idx++;
            } else // freq matching index must be determined if it is just higher than the freq at
                   // the previous point
            {
                start_indices_()(nextpointonrayindex, next_freq_idx)[0] =
                    currpointonrayindex; // stores point
                start_indices_()(nextpointonrayindex, next_freq_idx)[1] =
                    curr_freq_idx; // stores freq
                // for the sake of boundary conditions, these start_indices_ might be overwritten
                // later
                next_freq_idx++;
            }
        }
        // match remaining points as best as possible (with the outermost point)
        while (next_freq_idx <= max_freq_idx) {
            start_indices_()(nextpointonrayindex, next_freq_idx)[0] =
                currpointonrayindex;                                                // stores point
            start_indices_()(nextpointonrayindex, next_freq_idx)[1] = max_freq_idx; // stores freq
            next_freq_idx++;
        }
    }
}

///  Comoving solver setup helper function: Computes which lines are overlapping at the given point
///  (at the next point)
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] nextpoint: the point index of the next point
///   @param[in] is_upward_disc: whether or not to use the rightmost frequency indices as boundary
///   conditions
///   @note This step is necessary for determining the boundary conditions (we require to know which
///   lines overlap to determine those)
inline void Solver ::get_overlapping_lines(
    Model& model, const Size nextpoint, bool is_upward_disc) {
    // For checking whether the lines overlap, we simply check whether the frequency quadrature
    // overlaps Should be conceptually more simple than using the line centers and widths to figure
    // out whether they overlap
    Vector<unsigned char>& line_quad_discdir_overlap = line_quad_discdir_overlap_();
    for (Size lineidx = 0; lineidx < model.parameters->nlines(); lineidx++) {
        line_quad_discdir_overlap[lineidx] = false; // initialize to no overlap
    }

    Vector<Size>& line_count =
        line_count_(); // counts the number of quadratures encountered for each line

    if (is_upward_disc) {
        for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
            // get corresponding line of the sorted frequency index
            const Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(
                nextpoint, freqidx); // arbitrary frame, as within a single point, one does not need
                                     // to care about doppler shifts
            const Size lineidx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            // increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            // if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx] == model.parameters->nquads()) {
                line_count[lineidx] =
                    0; // reset the index (as counting this line is no longer needed)
                // Then if no other lines are currently being counted (i.e. entire count vectors is
                // zero),
                //  we conclude that the current continguous range has ended
                //  TODO SIMPLIFY: LOGIC USING SET AND TOTAL COUNT
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(),
                        [](int i) { return i == 0; })) {
                    // just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx] = true;
                }
            }
        }
    } else { // long for more easily handling reverse loops
        for (long freqidx = model.parameters->nfreqs() - 1; freqidx >= 0; freqidx--) {
            // get corresponding line of the sorted frequency index
            const Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(
                nextpoint, freqidx); // arbitrary frame, as within a single point, one does not need
                                     // to care about doppler shifts
            const Size lineidx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            // increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            // if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx] == model.parameters->nquads()) {
                line_count[lineidx] =
                    0; // reset the index (as counting this line is no longer needed)
                // Then if no other lines are currently being counted (i.e. entire count vectors is
                // zero),
                //  we conclude that the current continguous range has ended
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(),
                        [](int i) { return i == 0; })) {
                    // just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx] = true;
                }
            }
        }
    }
    // TODO: check line_count to be uniformly 0!! We did not reinitialize this variable...
}

///  Comoving solver setup helper function: Computes the continguous ranges spanned by the line
///  quadratures (at the current point)
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] curr_point: the point index of the current point
///   @param[in] is_upward_disc: whether or not to use the rightmost frequency indices as boundary
///   conditions
///   @param[in] curr_shift: doppler shift at the current point
///   @warning
///  @note: In this version, we exclude the 2 farthest frequency points due to boundary conditions
///  @note: Requires the shift, as we will use these computed ranges to determine overlap with
///  boundary frequencies
inline void Solver ::get_line_ranges(
    Model& model, const Size curr_point, bool is_upward_disc, Real curr_shift) {
    Vector<Real>& left_bound  = left_bound_();  // Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound = right_bound_(); // Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges           = nb_ranges_();   // contains the number of ranges
    Vector<Size>& left_bound_index = left_bound_index_(); // Specifies the corresponding freq index
                                                          // to the left bounds of the ranges
    Vector<Size>& right_bound_index =
        right_bound_index_(); // Specifies the corresponding freq index to the right bounds of the
                              // ranges

    nb_ranges = 0; // by default, set to 0.
    // First, we specify which quadratues are counted (every one except the outmost 2 ones depending
    // on discretization direction)
    Vector<Size>& quad_range_weight = quad_range_weight_();
    Size& tot_quad_range_weight     = tot_quad_range_weight_();
    tot_quad_range_weight =
        model.parameters->nquads() - 2; // Number of used quadratures is the same for both cases, so
                                        // we do not need to worry about this

    // TODO: maybe replace this implementation with the following:
    // const Size tot_quad_range_weight = model.parameters->nquads()-2;
    // Size tot_quad_count = 0;
    // std::set<Size> encountered_lines = std::set<Size>();
    //... In this, the set might make it a bit more efficient ...

    // We exclude the farthest 2 frequency points for each line, as these correspond to boundary
    // conditions
    if (is_upward_disc) { // then everything except the final two quads may be used for determining
                          // the range
        Size quadidx = 0;
        while (quadidx < model.parameters->nquads() - 2) {
            quad_range_weight[quadidx] = 1;
            quadidx++;
        }

        while (quadidx < model.parameters->nquads()) {
            quad_range_weight[quadidx] = 0;
            quadidx++;
        }
    } else { // then everything except the first two quads may be used for determining the range
        Size quadidx = 0;
        while (quadidx < 2) {
            quad_range_weight[quadidx] = 0;
            quadidx++;
        }

        while (quadidx < model.parameters->nquads()) {
            quad_range_weight[quadidx] = 1;
            quadidx++;
        }
    }

    // Dummy initialization, as we do not know a priori what the first counted quadrature will be
    Real leftbound     = 0.0;
    Real rightbound    = 0.0;
    Size leftboundidx  = 0;
    Size rightboundidx = 0;

    bool leftbound_specified =
        false; // Specifies whether the left bound for the current range has already been set

    Vector<Size>& line_count =
        line_count_(); // counts the number of quadratures encountered for each line
    // ASSUMED AT THIS POINT TO BE INITIALIZED TO 0.

    for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
        // get corresponding line of the sorted frequency index
        const Size unsorted_freqidx =
            model.radiation.frequencies.corresponding_nu_index(curr_point, freqidx);
        const Size lineidx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        const Size quadidx = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
        const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, freqidx)
                             * curr_shift; // in static frame

        if (!leftbound_specified
            && (quad_range_weight[quadidx])) { // specify left bound if not yet done so (and the
                                               // point is not a boundary point)
            leftbound_specified = true;
            leftbound           = curr_freq;
            leftboundidx        = freqidx;
        }
        rightbound    = curr_freq;
        rightboundidx = freqidx;

        // TODO: implementation using a set
        // //better implementation of checking whether we have counted all useful quadrature points
        // of a (maybe overlapping) line region. encountered_lines.insert(lineidx);
        // tot_quad_count+=quad_range_weight[quadidx];//the boundary points need not be counted,
        // //err, I accidentally forgot that we do not increment the tot_quad_count when approaching
        // bdy points (adding nonspecified bounds makes no sense) if
        // ((leftbound_specified)&&((encountered_lines.size()*tot_quad_range_weight) ==
        // tot_quad_count))
        // {
        //     //add range
        //     left_bound[nb_ranges]=leftbound;
        //     right_bound[nb_ranges]=rightbound;
        //     left_bound_index[nb_ranges]=leftboundidx;
        //     right_bound_index[nb_ranges]=rightboundidx;
        //     nb_ranges++;
        //
        //     //set flag for new left bound
        //     leftbound_specified=false;
        // }

        // increment number of freqs encountered if a non-boundary point is encountered (implicit
        // conversion from bool to int)
        line_count[lineidx] += quad_range_weight[quadidx];

        // if we have counted all relevant quads belonging to a specific line:
        if (line_count[lineidx] == tot_quad_range_weight) {
            line_count[lineidx] = 0;
            // Then if no other lines are currently being counted (i.e. entire count is zero),
            //  we conclude that the current continguous range has ended
            // TODO: make this thing scalable, using a set for the lines encoutered, and some
            // counter for total counts!
            if (std::all_of(
                    line_count.vec.begin(), line_count.vec.end(), [](int i) { return i == 0; })) {
                // add range
                left_bound[nb_ranges]        = leftbound;
                right_bound[nb_ranges]       = rightbound;
                left_bound_index[nb_ranges]  = leftboundidx;
                right_bound_index[nb_ranges] = rightboundidx;
                nb_ranges++;

                // set flag for new left bound
                leftbound_specified = false;
            }
        }
    }
    // TODO: check line_count to be uniformly 0!! As we do not reinitialize it anywhere (and we
    // assume it to be 0 when entering this and another function)
}

///  Comoving solver setup helper function: Sets all implicit boundary conditions, using the
///  computed frequency ranges for currpoint and overlaps for nextpoint. Uses the line ranges
///  computed by get_line_ranges to figure out which frequencies (at the next point) lie out of
///  bounds
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] nextpoint: the point index of the next point
///   @param[in] nextpointonrayindex: the ray position index of the next point
///   @param[in] shift_next: doppler shift at the next point
///   @param[in,out] multimap_freq_to_bdy_index: map which stores which boundary conditions have yet
///   to be matched
///   @param[in] is_upward_disc: whether or not to use the rightmost frequency indices as boundary
///   conditions
///   @warning Both Solver::get_overlapping_lines and Solver::get_line_ranges should be called
///   before calling this function.
inline void Solver ::set_implicit_boundary_frequencies(Model& model, const Size nextpoint,
    const Size nextpointonrayindex, Real shift_next,
    std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index, bool is_upward_disc) {
    Size curr_range_index                    = 0;
    Vector<Real>& left_bound                 = left_bound_();
    Vector<Real>& right_bound                = right_bound_();
    Vector<unsigned char>& line_quad_overlap = line_quad_discdir_overlap_();
    Size nb_ranges                           = nb_ranges_();
    Real left_bound_freq                     = left_bound[curr_range_index];
    Real right_bound_freq                    = right_bound[curr_range_index];

    if (is_upward_disc) {
        // The only difference between upward and downward direction lies in the location of the
        // forced boundary conditions
        for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
            Real next_freq = model.radiation.frequencies.sorted_nu(nextpoint, freqidx)
                           * shift_next; // in static frame
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);
            const Size lineidx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            const Size quadidx =
                model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];

            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx >= model.parameters->nquads() - 2 && !line_quad_overlap[lineidx]) {
                // Thus we add this frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(
                    next_freq, std::make_tuple(nextpointonrayindex, freqidx));
                continue;
            }
            // ranges are sorted, just like the freqs; therefore we can relatively efficently
            // compare the non-overlap
            //  by just advancing the index for the ranges
            while ((right_bound_freq < next_freq) && (curr_range_index < nb_ranges - 1)) {
                curr_range_index++;
                left_bound_freq  = left_bound[curr_range_index];
                right_bound_freq = right_bound[curr_range_index];
            }
            // The frequency lies outside the range if smaller than left_bound OR larger than
            // right_bound
            if (next_freq < left_bound_freq || next_freq > right_bound_freq) {
                // Thus we add this frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(
                    next_freq, std::make_tuple(nextpointonrayindex, freqidx));
            }
        }
    } else {
        for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
            Real next_freq = model.radiation.frequencies.sorted_nu(nextpoint, freqidx)
                           * shift_next; // in static frame
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);
            const Size lineidx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            const Size quadidx =
                model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx < 2 && !line_quad_overlap[lineidx]) {
                // Thus we add the frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(
                    next_freq, std::make_tuple(nextpointonrayindex, freqidx));
                continue;
            }
            // ranges are sorted, just like the freqs; therefore we can relatively efficently
            // compare the non-overlap
            while ((right_bound_freq < next_freq) && (curr_range_index < nb_ranges - 1)) {
                curr_range_index++;
                left_bound_freq  = left_bound[curr_range_index];
                right_bound_freq = right_bound[curr_range_index];
            }
            // The frequency lies outside the range if smaller than left_bound OR larger than
            // right_bound
            if (next_freq < left_bound_freq || next_freq > right_bound_freq) {
                // Thus we add the frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(
                    next_freq, std::make_tuple(nextpointonrayindex, freqidx));
            }
        }
    }
}

///  Comoving solver setup helper function: Matches the overlapping boundary conditions, using the
///  currently computed ranges.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] currpoint: the point index of the current point
///   @param[in] curr_point_on_ray_index: the ray position index of the current point
///   @param[in] curr_shift: doppler shift at the current point
///   @param[in,out] multimap_freq_to_bdy_index: map which stores which boundary conditions have yet
///   to be matched
///   @warning Should be called after Solver::set_implicit_boundary_frequencies
inline void Solver ::match_overlapping_boundary_conditions(Model& model, const Size currpoint,
    const Size curr_point_on_ray_index, const Real curr_shift,
    std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index) {
    Vector<Real>& left_bound  = left_bound_();  // Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound = right_bound_(); // Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges           = nb_ranges_();   // contains the number of ranges
    Vector<Size>& left_bound_index = left_bound_index_(); // Specifies the corresponding freq index
                                                          // to the left bounds of the ranges
    Vector<Size>& right_bound_index =
        right_bound_index_(); // Specifies the corresponding freq index to the right bounds of the
                              // ranges

    Size curr_range_index = 0;
    Size curr_freq_index =
        left_bound_index[curr_range_index]; // contains the current frequency index
    // Assumption: at least a single line exists -> a single range exists
    Real curr_range_min  = left_bound[curr_range_index];
    Real curr_range_max  = right_bound[curr_range_index];
    Real curr_range_freq = curr_range_min;

    // By just using the builtin lower_bound, one can easily find the boundary conditions
    // overlapping with the current ranges
    std::multimap<Real, std::tuple<Size, Size>>::iterator it =
        multimap_freq_to_bdy_index.lower_bound(curr_range_min);
    // Using a while loop, as we simulataneously iterate over the boundary conditions and over the
    // line ranges
    while (it != multimap_freq_to_bdy_index.end() && curr_range_index < nb_ranges) {
        const Real bdy_freq = it->first; // get boundary freq
        const Size bdy_point_on_ray_idx =
            std::get<0>(it->second);                       // get corresponding point index on ray
        const Size bdy_freq_idx = std::get<1>(it->second); // get corresponding freq index
        // check if bdy freq is smaller than the right bound
        if (bdy_freq > curr_range_max) {
            // The boundary frequency lies further than the current line, so use the next line
            // instead
            curr_range_index++;
            if (curr_range_index == nb_ranges) { // no next line exist, so all remaining boundary
                                                 // frequencies lie outside the last range
                return;
            }
            // update info about ranges
            curr_range_min  = left_bound[curr_range_index];
            curr_range_max  = right_bound[curr_range_index];
            curr_range_freq = curr_range_min;
            curr_freq_index =
                left_bound_index[curr_range_index]; // stores the freq index of curr point
            it = multimap_freq_to_bdy_index.lower_bound(
                curr_range_min); // and adjust iterator accordingly
            continue;
        }
        // the bdy frequency now lies within the bounds, so just linearly iterate until we find the
        // correct freqs of the ranges (static frame) to interpolate with

        // Find which freq bounds at curr_point exactly correspond to the bdy freq (enclosing it
        // from the right)
        while (bdy_freq > curr_range_freq && curr_range_freq < curr_range_max) {
            curr_freq_index++;
            curr_range_freq = model.radiation.frequencies.sorted_nu(currpoint, curr_freq_index)
                            * curr_shift; // in static frame
        }

        // In this section: I have commented a few different options for interpolating the
        // boundary conditions

        /// zeroth order interpolation for the boundary frequency
        // Size left_curr_freq_idx=curr_freq_index-1;
        // Size right_curr_freq_idx=curr_freq_index;
        //
        // //the curr range freq should now be larger/equal to the bdy freq
        // //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as
        // index for interpolation, so use curr_index as left bound instead if
        // (curr_freq_index==left_bound_index[curr_range_index])
        // {
        //   left_curr_freq_idx=curr_freq_index;
        //   right_curr_freq_idx=curr_freq_index+1;
        // }
        //
        // const Real left_freq = model.radiation.frequencies.sorted_nu(currpoint,
        // left_curr_freq_idx)*curr_shift; const Real right_freq =
        // model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)*curr_shift; Size
        // zeroth_order_freq = 0;
        //
        // //get closest freq
        // if (bdy_freq-left_freq<0.5*(right_freq-left_freq))
        // {
        //     zeroth_order_freq = left_curr_freq_idx;
        // }
        // else
        // {
        //     zeroth_order_freq = right_curr_freq_idx;
        // }

        /// First order interpolation for the boundary frequency
        Size left_curr_freq_idx  = curr_freq_index - 1;
        Size right_curr_freq_idx = curr_freq_index;

        // the curr range freq should now be larger/equal to the bdy freq
        // In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as
        // index for interpolation, so use curr_index as left bound instead
        if (curr_freq_index == left_bound_index[curr_range_index]) {
            left_curr_freq_idx  = curr_freq_index;
            right_curr_freq_idx = curr_freq_index + 1;
        }

        const Real deltafreq =
            (model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)
                - model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx))
            * curr_shift; // in static frame

        /// Second order interpolation
        // Size left_curr_freq_idx=curr_freq_index-2;
        // Size middle_curr_freq_idx=curr_freq_index-1;
        // Size right_curr_freq_idx=curr_freq_index;
        // //the curr range freq should now be larger/equal to the bdy freq
        // //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as
        // index for interpolation, so use curr_index as left bound instead if
        // (curr_freq_index==left_bound_index[curr_range_index]||middle_curr_freq_idx==left_bound_index[curr_range_index])
        // {
        //     left_curr_freq_idx=left_bound_index[curr_range_index];
        //     middle_curr_freq_idx=left_bound_index[curr_range_index]+1;
        //     right_curr_freq_idx=left_bound_index[curr_range_index]+2;
        // }
        //
        // //first compute coefficents for the second order accurate freq derivative for the
        // explicit part const Real dfreqsmall=(model.radiation.frequencies.sorted_nu(currpoint,
        // right_curr_freq_idx)-model.radiation.frequencies.sorted_nu(currpoint,
        // middle_curr_freq_idx))
        //             *curr_shift;//in static frame
        // const Real dfreqlarge=(model.radiation.frequencies.sorted_nu(currpoint,
        // right_curr_freq_idx)-model.radiation.frequencies.sorted_nu(currpoint,
        // left_curr_freq_idx))
        //             *curr_shift;//in static frame

        /// Setting the interpolation coefficients for the boundary condition

        // Set frequency derivative correctly

        // Note: the actual factor with which to multiply is Δτ^2/(1-exp(-Δτ)-Δτ*exp(-Δτ)), which in
        // the limit Δτ→0 corresponds to 2 (wait, I did forget a factor exp(-Δτ) in front due to
        // optical depth itself) So in total, the multiplication factor is
        // exp(-Δτ)*Δτ^2/(1-exp(-Δτ)-Δτ*exp(-Δτ)) Set explicit coefficents to TWICE the normal 1/Δν
        // value (as we are only treating the explicit part)

        /// zeroth order interpolation
        // //set starting indices correctly for zeroth order interpolation
        // start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[0]=curr_point_on_ray_index;
        // start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1]=zeroth_order_freq;
        //
        // dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // //as the zeroth order does not use the freq derivative, the indices here do not matter
        // dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;

        // Set starting index correctly (for first, second order interpolation)
        start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[0] = curr_point_on_ray_index;
        start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1] = left_curr_freq_idx;
        /// first order interpolation
        dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = -2.0 / deltafreq;
        dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 2.0 / deltafreq;
        dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;

        dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = left_curr_freq_idx;
        dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = right_curr_freq_idx;
        dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = right_curr_freq_idx;
        // as the coefficient for the last on is 0, no effort should be made to renumber that last
        // index

        /// second order interpolation (might be worse due to trying to fit some discontinouous
        /// thing with a higher order polynomial)
        // dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)
        // =-2.0*dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
        // dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)
        // =2.0*dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
        // dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)
        // =-dIdnu_coef3_curr_()(bdy_point_on_ray_idx,
        // bdy_freq_idx)-dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx);//curr point itself
        //
        // dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=middle_curr_freq_idx;
        // dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;

        // Implicit part coefs are set to 0. (as long as the corresponding indicies are in bounds, i
        // do not particularly care about the exact values)
        dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        // exact indices do not matter, as long as they are in bounds
        dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;
        dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;
        dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;

        // set dtau to approx 0
        delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx) =
            model.parameters
                ->comoving_min_dtau; // should be small enough, but not small enough to
                                     // crash my solver (due to factor /Δτ^2 in computation)
        // The source function is set to 0, as the optical depth increment Δτ≃0 anyway.
        S_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        S_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        // TODO: In case of large doppler shifts, we might need to rethink setting dtau→0; as this
        // will underestimate the optical depth (and influence of source) Concretely, we would need
        // to use the optical depth of the last segment, and somehow incorporate the interpolated
        // intensity of way earlier (and set the source function correctly)

        // pop value from iterator, as we have succesfully set this boundary condition
        it = multimap_freq_to_bdy_index.erase(it);
    }
}

///  Comoving solver setup helper function: for all remaining boundary points, this sets the
///  boundary condition correctly. Uses the start of the ray for determining the boundary
///  condition
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] initial_bdy: the point index of the of the start of the ray
///   @param[in] curr_shift: doppler shift at the start of the ray
///   @param[in,out] multimap_freq_to_bdy_index: map which stores which boundary conditions have yet
///   to be matched
inline void Solver ::set_initial_boundary_conditions(Model& model, const Size initial_bdy,
    const Real curr_shift,
    std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index) {
    std::multimap<Real, std::tuple<Size, Size>>::iterator it = multimap_freq_to_bdy_index.begin();
    while (it != multimap_freq_to_bdy_index.end()) {
        const Real bdy_freq = it->first; // get boundary freq (in static frame)
        const Size bdy_point_on_ray_idx =
            std::get<0>(it->second);                       // get corresponding point index on ray
        const Size bdy_freq_idx = std::get<1>(it->second); // get corresponding freq index

        // For using somewhat reasonable frequency indices
        Size curr_freq_idx = start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1];

        // Compute boundary intensity
        Real bdy_intensity = boundary_intensity(model, initial_bdy, bdy_freq);

        // Set the frequency derivative coefficients to 0
        dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;

        dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = curr_freq_idx;
        dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = curr_freq_idx;
        dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = curr_freq_idx;

        dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;
        dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = 0.0;

        dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;
        dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;
        dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_freq_idx;

        // Set S and dtau such that the resulting intensity will be (approximately) equal to the
        // initial boundary intensity
        delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx) =
            50.0; // should be large enough for 1-e^(-dtau) ≃ 1
        S_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_intensity;
        S_next_()(bdy_point_on_ray_idx, bdy_freq_idx) = bdy_intensity;

        // TODO: Hmm, if we were to compute S, Δτ (only using the last part), we might be able to
        // get better results in case of extreme velocity gradients (and limited quadrature
        // points). However, for this to actually have an effect on the computed results, the
        // densities must be high enough. In that regime, the model is just poorly sampled anyway,
        // so this might not be worth the time to implement.

        // pop value from iterator, as this boundary condition is now accounted for
        it = multimap_freq_to_bdy_index.erase(it);
    }
}

///  Comoving solver imager post-processing: interpolates the intensities at the requested
///  frequencies using the computed intensities at currpoint.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] image: the image to save the intensities to
///   @param[in] pixidx: pixel index of image to save intensity to
///   @param[in] currpoint: the point index of the current point
///   @param[in] curr_point_on_ray_index: the ray position index of the current point
///   @param[in] curr_shift: doppler shift of the current point
///   @param[in,out] multimap_freq_to_bdy_index: map which stores which frequencies have yet to be
///   interpolated, and to what image frequency index it corresponds to
inline void Solver ::interpolate_computed_comoving_intensities(Model& model, Image& image,
    Size pixidx, const Size currpoint, const Size curr_point_on_ray_index, const Real curr_shift,
    std::multimap<Real, Size>& multimap_image_freq_to_index) {
    Vector<Real>& left_bound  = left_bound_();  // Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound = right_bound_(); // Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges           = nb_ranges_();   // contains the number of ranges
    Vector<Size>& left_bound_index = left_bound_index_(); // Specifies the corresponding freq index
                                                          // to the left bounds of the ranges
    Vector<Size>& right_bound_index =
        right_bound_index_(); // Specifies the corresponding freq index to the right bounds of the
                              // ranges

    Matrix<Real>& intensities = intensities_(); // contains the computed intensities at each point

    Size curr_range_index = 0;
    Size curr_freq_index =
        left_bound_index[curr_range_index]; // contains the current frequency index
    // Assumption: at least a single line exists -> a single range exists
    Real curr_range_min  = left_bound[curr_range_index];
    Real curr_range_max  = right_bound[curr_range_index];
    Real curr_range_freq = curr_range_min;

    // By just using the builtin lower_bound, one can easily find the boundary conditions
    // overlapping with the current ranges
    std::multimap<Real, Size>::iterator it =
        multimap_image_freq_to_index.lower_bound(curr_range_min);
    // Using a while loop, as we simulataneously iterate over the boundary conditions and over the
    // line ranges
    while (it != multimap_image_freq_to_index.end() && curr_range_index < nb_ranges) {
        const Real bdy_freq     = it->first;  // get image freq (already in static frame)
        const Size bdy_freq_idx = it->second; // get corresponding freq index for image

        // check if bdy freq is smaller than the right bound
        if (bdy_freq > curr_range_max) {
            // The boundary frequency lies further than the current line, so use the next line
            // instead
            curr_range_index++;
            if (curr_range_index == nb_ranges) { // no next line exist, so all remaining boundary
                                                 // frequencies lie outside the last range
                return;
            }
            // update info about ranges
            curr_range_min  = left_bound[curr_range_index];
            curr_range_max  = right_bound[curr_range_index];
            curr_range_freq = curr_range_min;
            curr_freq_index =
                left_bound_index[curr_range_index]; // stores the freq index of curr point
            it = multimap_image_freq_to_index.lower_bound(
                curr_range_min); // and adjust iterator accordingly
            continue;
        }
        // the bdy frequency now lies within the bounds, so just linearly iterate until we find the
        // correct freqs of the ranges (static frame) to interpolate with

        // Find which freq bounds at curr_point exactly correspond to the bdy freq (enclosing it
        // from the right)
        while (bdy_freq > curr_range_freq && curr_range_freq < curr_range_max) {
            curr_freq_index++;
            curr_range_freq = model.radiation.frequencies.sorted_nu(currpoint, curr_freq_index)
                            * curr_shift; // in static frame
        }

        /// First order interpolation for the boundary frequency
        Size left_curr_freq_idx  = curr_freq_index - 1;
        Size right_curr_freq_idx = curr_freq_index;

        // the curr range freq should now be larger/equal to the bdy freq
        // In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as
        // index for interpolation, so use curr_index as left bound instead
        if (curr_freq_index == left_bound_index[curr_range_index]) {
            left_curr_freq_idx  = curr_freq_index;
            right_curr_freq_idx = curr_freq_index + 1;
        }

        const Real deltafreq =
            (model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)
                - model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx))
            * curr_shift; // in static frame

        const Real delta_I = intensities(curr_point_on_ray_index, right_curr_freq_idx)
                           - intensities(curr_point_on_ray_index, left_curr_freq_idx);

        const Real interpolated_intensity =
            intensities(curr_point_on_ray_index, left_curr_freq_idx)
            + delta_I / deltafreq
                  * (bdy_freq
                      - model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx)
                            * curr_shift);

        image.I(pixidx, bdy_freq_idx) = interpolated_intensity;

        // pop value from iterator, as we have succesfully interpolated the intensity for this
        // frequency
        it = multimap_image_freq_to_index.erase(it);
    }
}

///  Comoving solver setup helper function: Does the required setup for the comoving solver in the
///  forward ray direction. Also correctly incorporates the boundary conditions. This function
///  assumes we have already traced the ray and this maps all data required for the actual
///  computation
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] last_interesting_rayposidx: last ray position index needed for the comoving
///   computation
template <ApproximationType approx>
inline void Solver ::comoving_ray_bdy_setup_forward(Model& model, Size last_interesting_rayposidx) {
    // This assumes we use a forward ray, not a backward ray
    //  The main practical differences lie in the shift direction and the traversal direction of the
    //  ray
    Matrix<Real>& intensities = intensities_();

    Size first_index = first_(); // index of first point on ray
    Real shift_first = 2.0 - shift_()[first_()];
    // As we only need to compute stuff on the ray until the last_interesting_rayposidx, we will
    // start from there (instead of the last_() ray position index)
    Size rayposidx =
        last_interesting_rayposidx; // ray position index -> point index through nr_()[rayposidx]
    // Initialize intensities at the start of the ray
    for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
        intensities(first_index, freqid) = boundary_intensity(model, nr_()[first_index],
            model.radiation.frequencies.sorted_nu(nr_()[first_index], freqid) * shift_first);
    }

    // dummy initialization; just need some references to Reals for get_eta_and_chi
    Real eta_next = 0.0;
    Real chi_next = 0.0;
    Real eta_curr = 0.0;
    Real chi_curr = 0.0;

    // helper values for computing the freq derivative terms
    Real dfreqsmall = 0.0;
    Real dfreqlarge = 0.0;

    std::multimap<Real, std::tuple<Size, Size>>
        multimap_freq_to_bdy_index; // stores the boundary conditions
    // boundary conditions can overlap, so multimap it is

    bool computing_curr_opacity = true;

    // For determining the boundary conditions, we are going backwards; this seems strange, but this
    // results in an easier treatment of the spatial boundary conditions at the end
    while (rayposidx > first_()) {
        const Size nextpoint = nr_()[rayposidx];
        const Size currpoint = nr_()[rayposidx - 1];
        const Real dZ        = dZ_()[rayposidx - 1];

        const Real shift_next =
            2.0 - shift_()[rayposidx]; // shift is by default defined backwards for the forward ray
                                       // direction; so under the non-relativistic approx, we can
                                       // just do 2-shift to get the shift in the other direction
        const Real shift_curr =
            2.0
            - shift_()[rayposidx - 1]; // this is due to its meaning; it is the transformation from
                                       // the static to the comoving frame of that specific point

        // For the sake of accurately putting boundary conditions, this should point in the correct
        // direction; otherwise we might accidentally put boundary conditions somewhat closer to the
        // middle of a line
        const bool is_upward_disc = (shift_next >= shift_curr);

        // TODO: maybe refactor loop from this point onwards; should be the same for both forward
        // and backward rays

        // First, we match all frequency indices as good as possible. This is later used for
        // computing S, Δτ
        match_frequency_indices(model, nextpoint, currpoint, shift_next, shift_curr, rayposidx,
            rayposidx - 1, is_upward_disc); // O(nfreqs)

        bool compute_curr_opacity = true;
        // Now compute all default stuff, computing bogus for the non-matched indices (but as these
        // correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Real nextfreq = model.radiation.frequencies.sorted_nu(
                nextpoint, next_freq_idx); //=comoving frame freq at next point
            const Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq      = model.radiation.frequencies.sorted_nu(
                     currpoint, curr_freq_idx); //=comoving frame freq at curr point (different frame)
            // technically, we also need to read the point index from this start_indices_
            // (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to
            // currpoint at this moment in time
            //  Retrieving the line index from the sorted frequency index
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpoint, next_freq_idx);
            const Size nextlineidx =
                model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            // only useful for single line approx

            Real dtau, Snext, Scurr;
            Real chicurr, chinext; // dummy stuff

            // The new method for computing S, dtau can handle moving frequencies without issues
            // TODO: figure out a way to only selectively need to compute the opacities twice (for
            // each frequency seperately); i.e replace: compute_curr_opacity = true; with something
            // more elaborate
            compute_curr_opacity = true;
            compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq,
                nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext,
                Scurr, Snext);
            S_next_()(rayposidx, next_freq_idx) = Snext;
            S_curr_()(rayposidx, next_freq_idx) = Scurr;
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            dtau = std::max(dtau, model.parameters->comoving_min_dtau);
            delta_tau_()(rayposidx, next_freq_idx) = dtau;
        }

        // This part precomputes the frequency derivative coefficients (and corresponding indices),
        // conveniently ignoring the fact that for some outermost line quadrature frequencies, we
        // must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        // The discretization is evidently depends on freq discretization direction
        if (is_upward_disc) {
            // Do not forget to exclude the last two frequencies, as these cannot have any freq der
            // computation (should definitely be bdy conditions). NOTE: The solver cannot deal with
            // nfreqs < 3. When this happens, we only have boundary conditions.
            for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs() - 2;
                 next_freq_idx++) {
                Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
                // Modulo operator seems weird, but for the last few freqs, the corresponding freq
                // might be the last one of currpoint == nfreqs()-1 If we want to prevent oob
                // accesses and divides by zero, we need to remap the oob indices (my choice is just
                // using modulo). The exact remapped indices should not matter, as all of this will
                // be overwritten when determining the boundary indices
                Size curr_freq_idxp1 = (curr_freq_idx + 1) % model.parameters->nfreqs(); // index+1
                Size curr_freq_idxp2 = (curr_freq_idx + 2) % model.parameters->nfreqs(); // index+2
                // first compute coefficents for the second order accurate freq derivative for the
                // explicit part
                dfreqsmall =
                    (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp2)
                        - model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp1))
                    * shift_curr;
                dfreqlarge = (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp2)
                                 - model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx))
                           * shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_curr_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_curr_()(rayposidx, next_freq_idx); // at curr_freq_idx itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx) = curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx) = curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx) = curr_freq_idx;

                // And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall =
                    (model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx + 2)
                        - model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx + 1))
                    * shift_next;
                dfreqlarge = (model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx + 2)
                                 - model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx))
                           * shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_next_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_next_()(rayposidx, next_freq_idx); // at next_freq_idx itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx) = next_freq_idx + 2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx) = next_freq_idx + 1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx) = next_freq_idx;
            }
        } else {
            // Do not forget to exclude the first two frequencies, as these cannot have any freq der
            // computation (should definitely be bdy conditions). NOTE: The solver cannot deal with
            // nfreqs < 3. When this happens, we only have boundary conditions.
            for (Size next_freq_idx = 2; next_freq_idx < model.parameters->nfreqs();
                 next_freq_idx++) {
                Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
                // Modulo operator seems weird, but for the last few freqs, the corresponding freq
                // might be the last one of currpoint == nfreqs()-1 If we want to prevent oob
                // accesses and divides by zero, we need to remap the oob indices (my choice is just
                // using modelo) The exact remapped indices should not matter, as all of this will
                // be overwritten when determining the boundary indices

                // MODULUS FOR 'negative' (more like overflowing) unsigned integer is in general not
                // that easy to compute. However, index only goes ever so slightly negative
                // ('-1,-2'), thus (mod+i)%mod is sufficient for the computation (mod>=2)
                Size curr_freq_idxm1 = (curr_freq_idx - 1 + model.parameters->nfreqs())
                                     % model.parameters->nfreqs(); // index-1
                Size curr_freq_idxm2 = (curr_freq_idx - 2 + model.parameters->nfreqs())
                                     % model.parameters->nfreqs(); // index-2
                // first compute coefficents for the second order accurate freq derivative for the
                // explicit part Note: I use minus signs in front, as we are now computing the first
                // derivative using points on the other side
                dfreqsmall =
                    -(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm1)
                        - model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm2))
                    * shift_curr;
                dfreqlarge =
                    -(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx)
                        - model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm2))
                    * shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_curr_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_curr_()(rayposidx, next_freq_idx); // at the curr_freq_idx itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx) = curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx) = curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx) = curr_freq_idx;

                // And now do the same (but simpler; less index management) for the implicit part.
                // Note: minus signs in front, as we are now computing the first derivative using
                // points on the other side
                dfreqsmall =
                    -(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx - 1)
                        - model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx - 2))
                    * shift_next;
                dfreqlarge =
                    -(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx)
                        - model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx - 2))
                    * shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_next_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_next_()(rayposidx, next_freq_idx); // at the next_freq_idx itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx) = next_freq_idx - 2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx) = next_freq_idx - 1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx) = next_freq_idx;
            }
        }

        // now start classifying the new boundary points
        // For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc,
            shift_curr); // for the current point, we definitely need the exact ranges (fiddled with
                         // to ensure enough bdy conditions)

        // and compute the overlap between lines on next_point to make sure we do not accidentally
        // add extra boundary conditions where lines might overlap
        get_overlapping_lines(model, nextpoint,
            is_upward_disc); // technically, we should compute which lines overlap only in this part

        // Then, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be
        // treated as boundary points
        set_implicit_boundary_frequencies(
            model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);

        // Finally, we check what boundary conditions need to be evaluated using the intensities at
        // curr_point
        match_overlapping_boundary_conditions(
            model, currpoint, rayposidx - 1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx--;
    }
    // After going through all points, the remaining boundary frequencies have not matched any
    // frequencies near any line, so we will use initial boundary conditions computed at the first
    // point to put on the ray note: we use the first point on the ray for determining the boundary
    // intensity, as the boundary condition originates from that point
    set_initial_boundary_conditions(
        model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
}

///  Comoving solver setup helper function: Does the required setup for the comoving solver in the
///  backward ray direction. Also correctly incorporates the boundary conditions. This function
///  assumes we have already traced the ray and this maps all data required for the actual
///  computation
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] first_interesting_rayposidx: first ray position index needed for the comoving
///   computation
template <ApproximationType approx>
inline void Solver ::comoving_ray_bdy_setup_backward(
    Model& model, Size first_interesting_rayposidx) {
    Matrix<Real>& intensities = intensities_();

    Size first_index = last_();               // index of first point on ray
    Real shift_first = shift_()[first_index]; // shift on backward ray is reverse of shift of
                                              // forward ray (which itself is the reverse operation
                                              // of mapping from comoving to static frame)
    Size rayposidx =
        first_interesting_rayposidx; // first index for which we want to know the intensity, thus
                                     // all indices before can be ignored
    for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
        intensities(first_index, freqid) = boundary_intensity(model, nr_()[first_index],
            model.radiation.frequencies.sorted_nu(nr_()[first_index], freqid) * shift_first);
    }

    // dummy initialization; just need some references to Reals for get_eta_and_chi
    Real eta_next = 0.0;
    Real chi_next = 0.0;
    Real eta_curr = 0.0;
    Real chi_curr = 0.0;

    // helper values for computing the freq derivative terms
    Real dfreqsmall = 0.0;
    Real dfreqlarge = 0.0;

    bool computing_curr_opacity = true;
    std::multimap<Real, std::tuple<Size, Size>>
        multimap_freq_to_bdy_index; // store all boundary condition frequencies, together with
                                    // auxiliary information (originating point and frequency index)
    // boundary condition frequencies can overlap, so multimap it is

    // In order to more easily determine the boundary conditions, we trace the ray backwards
    // (compared to the direction). This allows for a simpler treatment of the boundary conditions.
    while (rayposidx < last_()) {
        const Size nextpoint = nr_()[rayposidx];
        const Size currpoint = nr_()[rayposidx + 1];
        const Real dZ        = dZ_()[rayposidx]; // dZ is stored somewhat finnicky, in locations
                                                 // [first_(),last_()[ (so not including last_())

        const Real shift_next =
            shift_()[rayposidx]; // shift is by default defined backwards for the forward ray
                                 // direction; so under the non-relativistic approx, we can just do
                                 // 2-shift to get the shift in the other direction
        const Real shift_curr = shift_()[rayposidx + 1];

        // For the sake of accurately putting boundary conditions, this should point in the correct
        // direction; otherwise we might accidentally put boundary conditions somewhat closer to the
        // middle of a line
        const bool is_upward_disc = (shift_next >= shift_curr);

        // TODO: maybe refactor loop from this point onwards; should be the same for both forward
        // and backward rays

        // First, we match all frequency indices as good as possible, in order to be able to compute
        // source functions and optical depths
        match_frequency_indices(model, nextpoint, currpoint, shift_next, shift_curr, rayposidx,
            rayposidx + 1, is_upward_disc); // O(nfreqs)

        bool compute_curr_opacity = true;
        // Now compute all default stuff, computing bogus for the non-matched indices (but as these
        // correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Real nextfreq = model.radiation.frequencies.sorted_nu(
                nextpoint, next_freq_idx); //=comoving frame freq
            const Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq =
                model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx); //
            // Get the line index from the unsorted frequency index
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpoint, next_freq_idx);
            const Size nextlineidx =
                model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            Real dtau, Snext, Scurr;
            Real chicurr, chinext; // dummy initializations

            // The new method for computing S, dtau can handle moving frequencies without issues
            // TODO: figure out a way to only selectively need to compute the opacities twice (for
            // each frequency seperately); i.e replace: compute_curr_opacity = true; with something
            // more elaborate
            compute_curr_opacity = true;
            compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq,
                nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext,
                Scurr, Snext);
            S_next_()(rayposidx, next_freq_idx) = Snext;
            S_curr_()(rayposidx, next_freq_idx) = Scurr;
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            dtau = std::max(dtau, model.parameters->comoving_min_dtau);
            delta_tau_()(rayposidx, next_freq_idx) = dtau;
        }

        // This part precomputes the frequency derivative coefficients (and corresponding indices),
        // conveniently ignoring the fact that for some outermost line quadrature frequencies, we
        // must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        // The discretization is evidently depends on freq disc direction
        if (is_upward_disc) {
            // Do not forget to exclude the last two frequencies, as these cannot have any freq der
            // computation (should definitely bdy conditions)
            for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs() - 2;
                 next_freq_idx++) {
                Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
                // Modulo operator seems weird, but for the last few freqs, the corresponding freq
                // might be the last one of currpoint == nfreqs()-1 If we want to prevent oob
                // accesses and divides by zero, we need to remap the oob indices (my choice is just
                // using modelo). The exact remapped indices should not matter, as all of this will
                // be overwritten when determining the boundary indices
                Size curr_freq_idxp1 = (curr_freq_idx + 1) % model.parameters->nfreqs(); // index+1
                Size curr_freq_idxp2 = (curr_freq_idx + 2) % model.parameters->nfreqs(); // index+2
                // first compute coefficents for the second order accurate freq derivative for the
                // explicit part
                dfreqsmall = (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)
                                 - model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
                           * shift_curr;
                dfreqlarge = (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)
                                 - model.radiation.frequencies.nu(currpoint, curr_freq_idx))
                           * shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_curr_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_curr_()(
                        rayposidx, next_freq_idx); // corresponds to the frequency index at the
                                                   // current point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx) = curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx) = curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx) = curr_freq_idx;

                // And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall = (model.radiation.frequencies.nu(nextpoint, next_freq_idx + 2)
                                 - model.radiation.frequencies.nu(nextpoint, next_freq_idx + 1))
                           * shift_next;
                dfreqlarge = (model.radiation.frequencies.nu(nextpoint, next_freq_idx + 2)
                                 - model.radiation.frequencies.nu(nextpoint, next_freq_idx))
                           * shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_next_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_next_()(
                        rayposidx, next_freq_idx); // corresponds to the frequency index at the next
                                                   // point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx) = next_freq_idx + 2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx) = next_freq_idx + 1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx) = next_freq_idx;
            }
        } else {
            // Do not forget to exclude the first two frequencies, as these cannot have any freq der
            // computation (should definitely bdy conditions)
            for (Size next_freq_idx = 2; next_freq_idx < model.parameters->nfreqs();
                 next_freq_idx++) {
                Size curr_freq_idx = start_indices_()(rayposidx, next_freq_idx)[1];
                // Modulo operator seems weird, but for the last few freqs, the corresponding freq
                // might be the last one of currpoint == nfreqs()-1 If we want to prevent oob
                // accesses and divides by zero, we need to remap the oob indices (my choice is just
                // using modelo) The exact remapped indices should not matter, as all of this will
                // be overwritten when determining the boundary indices

                // MODULUS FOR 'negative' (more like overflowing) unsigned integer is in general not
                // that easy to compute However, as i know it only goes ever so slightly negative
                // ('-1,-2'): (mod+i)%mod is sufficient for the computation (mod>=1,2)
                Size curr_freq_idxm1 = (curr_freq_idx - 1 + model.parameters->nfreqs())
                                     % model.parameters->nfreqs(); // index-1
                Size curr_freq_idxm2 = (curr_freq_idx - 2 + model.parameters->nfreqs())
                                     % model.parameters->nfreqs(); // index-2
                // first compute coefficents for the second order accurate freq derivative for the
                // explicit part Note: minus signs in front, as we are now computing the first
                // derivative using points on the other side (switched order of term, making sign
                // clearer; just compare to upward disc version)
                dfreqsmall = -(model.radiation.frequencies.nu(currpoint, curr_freq_idxm1)
                                 - model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                           * shift_curr;
                dfreqlarge = -(model.radiation.frequencies.nu(currpoint, curr_freq_idx)
                                 - model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                           * shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_curr_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_curr_()(rayposidx, next_freq_idx); // curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx) = curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx) = curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx) = curr_freq_idx;

                // And now do the same (but simpler; less index management) for the implicit part
                // Note: minus signs in front, as we are now computing the first derivative using
                // points on the other side
                dfreqsmall = -(model.radiation.frequencies.nu(nextpoint, next_freq_idx - 1)
                                 - model.radiation.frequencies.nu(nextpoint, next_freq_idx - 2))
                           * shift_next;
                dfreqlarge = -(model.radiation.frequencies.nu(nextpoint, next_freq_idx)
                                 - model.radiation.frequencies.nu(nextpoint, next_freq_idx - 2))
                           * shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =
                    -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =
                    dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =
                    -dIdnu_coef3_next_()(rayposidx, next_freq_idx)
                    - dIdnu_coef2_next_()(rayposidx, next_freq_idx); // curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx) = next_freq_idx - 2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx) = next_freq_idx - 1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx) = next_freq_idx;
            }
        }

        // now start classifying the new boundary points
        // For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc,
            shift_curr); // for the current point, we definitely need the exact ranges (fiddled with
                         // to ensure enough bdy conditions)

        // and compute the overlap between lines on next_point
        get_overlapping_lines(model, nextpoint,
            is_upward_disc); // technically, we should compute which lines overlap only in this part

        // Now, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be
        // treated as boundary points
        set_implicit_boundary_frequencies(
            model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);

        // And we check what boundary conditions need to be evaluated using the intensities at
        // curr_point
        match_overlapping_boundary_conditions(
            model, currpoint, rayposidx + 1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx++;
    }
    // After going through all points, the remaining boundary frequencies have not matched any
    // frequencies near any line, so we will use initial boundary conditions computed at the first
    // point to put on the ray note: nr_()[first_index] should be point on boundary of domain; in
    // this way, we can simply add the corresponding boundary conditions
    set_initial_boundary_conditions(
        model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
}

///  Solver for the non-approximate comoving shortcharacteristics method. This is a first order in
///  space, second order in frequency accurate solver
///   @param[in] model: the model to apply the comoving solver to
///   @note This solver might suffer a bit from spatial inaccuracies, as the traced rays might not
///   go exactly through the positions
template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_comoving_order_2_sparse(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        accelerated_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k) = 0.0;
            }
        })
    }

    std::cout << "hnrays: " << model.parameters->hnrays() << std::endl;
    // For each ray, solve the radiative transfer equation in a comoving manner.

    std::vector<Size> cum_n_points_to_trace_ray_through(model.parameters->hnrays() + 1, 0);
    Size counter = 0; // will contain the total number of points to trace through
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
        counter += points_to_trace_ray_through[rr].size();
        cum_n_points_to_trace_ray_through[rr + 1] = counter;
    }

    // Parallelization over both rays and ray directions, as otherwise the load balancing might be a
    // bit iffy for smaller models
    accelerated_for(totalidx, counter, {
        auto cum_n_points_pointer = std::upper_bound(cum_n_points_to_trace_ray_through.begin(),
            cum_n_points_to_trace_ray_through.end(), totalidx);
        Size rr =
            std::distance(cum_n_points_to_trace_ray_through.begin(), cum_n_points_pointer) - 1;
        Size rayidx             = totalidx - cum_n_points_to_trace_ray_through[rr];
        const Size o            = points_to_trace_ray_through[rr][rayidx];
        const double dshift_max = get_dshift_max(model, o);
        solve_comoving_order_2_sparse<approx, use_adaptive_directions>(model, o, rr, rayidx, dshift_max);
    })
}

///  Solver for the non-approximate comoving methods. Computes the intensities for a single step on
///  the ray. Assumes the datastructures to be setup correctly.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] rayposidx: ray position inddex of the next point
// CHANGED ///   @param[in] rayidx: index of ray to trace in the given direction
///   @param[in] o: position index of ray origin used to trace the ray
///   @param[in] rr: ray direction index to check whether to save the current intensities (if the
///   ray lies closest in the current direction to the point)
///   @param[in] is_upward_disc: determines the ordering of the implicit parts of the computation
///   @param[in] forward_ray: whether or not the ray is traversed in the forward direction. changes
///   the sign of the doppler shift
template <bool use_adaptive_directions>
inline void Solver ::solve_comoving_single_step(Model& model, const Size rayposidx,
    const Size o, const Size rr, const bool is_upward_disc, const bool forward_ray) {
    Vector<Size>& nr = nr_(); // stores the exact point indices
    Vector<double>& shift =
        shift_(); // stores the shifts versus the static frame; warning: for forward rays, this
                  // should be changed to 2.0-shift; for backward rays, it is correct. This is due
                  // to mapping from comoving frame to observer frame
    Matrix<Vector<Size>>& start_indices =
        start_indices_(); // for getting the previous indices associated with computing the
                          // intensity for the next point
    Matrix<Real>& delta_tau   = delta_tau_();   // stores the optical depths
    Matrix<Real>& S_curr      = S_curr_();      // source function of current point
    Matrix<Real>& S_next      = S_next_();      // source function of next point
    Matrix<Real>& intensities = intensities_(); // computed intensities

    // Frequency derivative coefficients
    Matrix<Real>& dIdnu_coef1_curr = dIdnu_coef1_curr_();
    Matrix<Real>& dIdnu_coef2_curr = dIdnu_coef2_curr_();
    Matrix<Real>& dIdnu_coef3_curr = dIdnu_coef3_curr_();

    Matrix<Size>& dIdnu_index1_curr = dIdnu_index1_curr_();
    Matrix<Size>& dIdnu_index2_curr = dIdnu_index2_curr_();
    Matrix<Size>& dIdnu_index3_curr = dIdnu_index3_curr_();

    Matrix<Real>& dIdnu_coef1_next = dIdnu_coef1_next_();
    Matrix<Real>& dIdnu_coef2_next = dIdnu_coef2_next_();
    Matrix<Real>& dIdnu_coef3_next = dIdnu_coef3_next_();

    Matrix<Size>& dIdnu_index1_next = dIdnu_index1_next_();
    Matrix<Size>& dIdnu_index2_next = dIdnu_index2_next_();
    Matrix<Size>& dIdnu_index3_next = dIdnu_index3_next_();

    // rayposidx represents the next position on the ray
    const Size nextpointidx = nr[rayposidx];
    const Real next_shift   = (forward_ray)
                                ? 2.0 - shift[rayposidx]
                                : shift[rayposidx]; // absolute shifts, versus static frame!

    // Do the explicit part
    for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
        // Get the point indices
        const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx)[0];
        const Size curr_point_idx          = nr[curr_point_on_ray_index];
        const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx)[1];
        const Real curr_shift =
            (forward_ray) ? 2.0 - shift[curr_point_on_ray_index]
                          : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

        const Real deltanu =
            model.radiation.frequencies.nu(nextpointidx, next_freq_idx) * next_shift
            - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
        const Real dtau = delta_tau(rayposidx, next_freq_idx);

        const Real expl_term = (-expm1(-dtau) - dtau * exp(-dtau)) / dtau;

        const Real expl_freq_der =
            dIdnu_coef1_curr(rayposidx, next_freq_idx)
                * intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))
            + dIdnu_coef2_curr(rayposidx, next_freq_idx)
                  * intensities(
                      curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))
            + dIdnu_coef3_curr(rayposidx, next_freq_idx)
                  * intensities(
                      curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx));

        intensities(rayposidx, next_freq_idx) =
            intensities(curr_point_on_ray_index, curr_freq_idx) * exp(-dtau)
            + expl_term * S_curr(rayposidx, next_freq_idx) // source term
            + expl_term * expl_freq_der * deltanu / dtau;  // frequency derivative term
    }

    // Implicit part ordering depends on the discretization direction
    if (is_upward_disc) {
        // In case of an upward discretization, the implicit part needs to start from the largest
        // frequencies (where we have extra boundary conditions on the end). We also do a reverse
        // loop over unsigned integers, thus the indices need to be incremented by 1.
        for (Size next_freq_idx = model.parameters->nfreqs(); next_freq_idx > 0; next_freq_idx--) {
            const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx - 1)[0];
            const Size curr_point_idx          = nr[curr_point_on_ray_index];
            const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx - 1)[1];
            const Real curr_shift =
                (forward_ray)
                    ? 2.0 - shift[curr_point_on_ray_index]
                    : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

            const Real deltanu =
                model.radiation.frequencies.nu(nextpointidx, next_freq_idx - 1) * next_shift
                - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
            const Real dtau = delta_tau(rayposidx, next_freq_idx - 1);

            const Real impl_term = (dtau + expm1(-dtau)) / dtau;
            // only the parts of the other points (to subtract/add)
            const Real impl_freq_der =
                dIdnu_coef2_next(rayposidx, next_freq_idx - 1)
                    * intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx - 1))
                + dIdnu_coef3_next(rayposidx, next_freq_idx - 1)
                      * intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx - 1));

            intensities(rayposidx, next_freq_idx - 1) =
                (intensities(rayposidx, next_freq_idx - 1)
                    + impl_term * S_next(rayposidx, next_freq_idx - 1) // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - dIdnu_coef1_next(rayposidx, next_freq_idx - 1) * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    } else {
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx)[0];
            const Size curr_point_idx          = nr[curr_point_on_ray_index];
            const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx)[1];
            const Real curr_shift =
                (forward_ray)
                    ? 2.0 - shift[curr_point_on_ray_index]
                    : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

            const Real deltanu =
                model.radiation.frequencies.nu(nextpointidx, next_freq_idx) * next_shift
                - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
            const Real dtau = delta_tau(rayposidx, next_freq_idx);

            const Real impl_term = (dtau + expm1(-dtau)) / dtau;
            // only the parts of the other points (to subtract/add)
            const Real impl_freq_der =
                dIdnu_coef2_next(rayposidx, next_freq_idx)
                    * intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx))
                + dIdnu_coef3_next(rayposidx, next_freq_idx)
                      * intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx));

            intensities(rayposidx, next_freq_idx) =
                (intensities(rayposidx, next_freq_idx)
                    + impl_term * S_next(rayposidx, next_freq_idx) // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - dIdnu_coef1_next(rayposidx, next_freq_idx) * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    }

    if (intensity_origin[nextpointidx].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nextpointidx);
        Size rcur  = std::get<1>(valid_rcur);
    // Finally increment J if this ray lies closest to the point in the raydirection rr
    // if (closest_ray(rr, nextpointidx) == rayidx) {
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nextpointidx, rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            // Get the details about the line to which the current frequency belongs
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpointidx, freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nextpointidx, k) +=
                lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid);

            // Computing the ALI element
            const Size curr_point_on_ray_index = start_indices_()(rayposidx, freqid)[0];
            const Size curr_point_idx          = nr[curr_point_on_ray_index];
            const Size curr_freq_idx           = start_indices_()(rayposidx, freqid)[1];
            const Real curr_shift =
                (forward_ray)
                    ? 2.0 - shift[curr_point_on_ray_index]
                    : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

            const Real deltanu =
                model.radiation.frequencies.sorted_nu(nextpointidx, freqid) * next_shift
                - model.radiation.frequencies.sorted_nu(curr_point_idx, curr_freq_idx) * curr_shift;

            const Real dtau = delta_tau(rayposidx, freqid);

            const Real constant =
                lspec.linedata.A[k] * lspec.quadrature.weights[z]
                * model.geometry.rays
                      .weight[rr]; // for integrating over both frequency (weighted with profile
                                   // function) and angle. We also need the einstein A coefficient

            const Real source_term =
                (dtau + expm1(-dtau)); // '/dtau' has already been applied to simplify the fraction
            const Real lambdaterm =
                constant * source_term
                / (dtau - dIdnu_coef1_next(rayposidx, freqid) * source_term * deltanu);

            lspec.lambda.add_element(nextpointidx, k, nextpointidx, lambdaterm);
        }
    }
}

///  Computes the intensities for imaging using the comoving solver for a single step. Assumes the
///  datastructures to be initialized correctly.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] rayposidx: the ray position index of the next point
///   @param[in] is_upward_disc: determines the ordering of the implicit parts of the computation
///   @param[in] forward_ray: whether or not the ray is traversed in the forward direction. changes
///   the sign of the doppler shift
inline void Solver::solve_comoving_image_single_step(
    Model& model, const Size rayposidx, const bool is_upward_disc, const bool forward_ray) {
    Vector<Size>& nr = nr_(); // stores the exact point indices
    Vector<double>& shift =
        shift_(); // stores the shifts versus the static frame; warning: for forward rays, this
                  // should be changed to 2.0-shift; for backward rays, it is correct
    Matrix<Vector<Size>>& start_indices =
        start_indices_(); // for getting the previous indices associated with computing the
                          // intensity for the next point
    Matrix<Real>& delta_tau   = delta_tau_();   // stores the optical depths
    Matrix<Real>& S_curr      = S_curr_();      // source function at the current point
    Matrix<Real>& S_next      = S_next_();      // source function at the next point
    Matrix<Real>& intensities = intensities_(); // computed intensities

    // Frequency derivative coefficients
    Matrix<Real>& dIdnu_coef1_curr = dIdnu_coef1_curr_();
    Matrix<Real>& dIdnu_coef2_curr = dIdnu_coef2_curr_();
    Matrix<Real>& dIdnu_coef3_curr = dIdnu_coef3_curr_();

    Matrix<Size>& dIdnu_index1_curr = dIdnu_index1_curr_();
    Matrix<Size>& dIdnu_index2_curr = dIdnu_index2_curr_();
    Matrix<Size>& dIdnu_index3_curr = dIdnu_index3_curr_();

    Matrix<Real>& dIdnu_coef1_next = dIdnu_coef1_next_();
    Matrix<Real>& dIdnu_coef2_next = dIdnu_coef2_next_();
    Matrix<Real>& dIdnu_coef3_next = dIdnu_coef3_next_();

    Matrix<Size>& dIdnu_index1_next = dIdnu_index1_next_();
    Matrix<Size>& dIdnu_index2_next = dIdnu_index2_next_();
    Matrix<Size>& dIdnu_index3_next = dIdnu_index3_next_();

    // rayposidx represents the next position on the ray
    const Size nextpointidx = nr[rayposidx];
    const Real next_shift   = (forward_ray)
                                ? 2.0 - shift[rayposidx]
                                : shift[rayposidx]; // absolute shifts, versus static frame!

    // Do the explicit part
    for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
        // Get the point indices
        const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx)[0];
        const Size curr_point_idx          = nr[curr_point_on_ray_index];
        const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx)[1];
        const Real curr_shift =
            (forward_ray) ? 2.0 - shift[curr_point_on_ray_index]
                          : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

        const Real deltanu =
            model.radiation.frequencies.nu(nextpointidx, next_freq_idx) * next_shift
            - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
        const Real dtau = delta_tau(rayposidx, next_freq_idx);

        const Real expl_term = (-expm1(-dtau) - dtau * exp(-dtau)) / dtau;

        const Real expl_freq_der =
            dIdnu_coef1_curr(rayposidx, next_freq_idx)
                * intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))
            + dIdnu_coef2_curr(rayposidx, next_freq_idx)
                  * intensities(
                      curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))
            + dIdnu_coef3_curr(rayposidx, next_freq_idx)
                  * intensities(
                      curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx));

        // Add every part together: past intensity, source term and frequency shift term
        intensities(rayposidx, next_freq_idx) =
            intensities(curr_point_on_ray_index, curr_freq_idx) * exp(-dtau)
            + expl_term * S_curr(rayposidx, next_freq_idx) // source term
            + expl_term * expl_freq_der * deltanu / dtau;  // frequency derivative term
    }

    // Implicit part ordering depends on the discretization direction
    if (is_upward_disc) {
        // In case of an upward discretization, the implicit part needs to start from the largest
        // frequencies (where we have extra boundary conditions on the end). We also do a reverse
        // loop over unsigned integers, thus the indices need to be incremented by 1.
        for (Size next_freq_idx = model.parameters->nfreqs(); next_freq_idx > 0; next_freq_idx--) {
            const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx - 1)[0];
            const Size curr_point_idx          = nr[curr_point_on_ray_index];
            const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx - 1)[1];
            const Real curr_shift =
                (forward_ray)
                    ? 2.0 - shift[curr_point_on_ray_index]
                    : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

            const Real deltanu =
                model.radiation.frequencies.nu(nextpointidx, next_freq_idx - 1) * next_shift
                - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
            const Real dtau = delta_tau(rayposidx, next_freq_idx - 1);

            const Real impl_term = (dtau + expm1(-dtau)) / dtau;
            const Real impl_freq_der =
                dIdnu_coef2_next(rayposidx, next_freq_idx - 1)
                    * intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx - 1))
                + dIdnu_coef3_next(rayposidx, next_freq_idx - 1)
                      * intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx - 1));

            // Do the implicit part of the computation: add current intensity and source term, then
            // divide by the (1-freq derivative) term
            intensities(rayposidx, next_freq_idx - 1) =
                (intensities(rayposidx, next_freq_idx - 1)
                    + impl_term * S_next(rayposidx, next_freq_idx - 1) // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - dIdnu_coef1_next(rayposidx, next_freq_idx - 1) * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    } else {
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Size curr_point_on_ray_index = start_indices_()(rayposidx, next_freq_idx)[0];
            const Size curr_point_idx          = nr[curr_point_on_ray_index];
            const Size curr_freq_idx           = start_indices_()(rayposidx, next_freq_idx)[1];
            const Real curr_shift =
                (forward_ray)
                    ? 2.0 - shift[curr_point_on_ray_index]
                    : shift[curr_point_on_ray_index]; // absolute shift, versus static frame

            const Real deltanu =
                model.radiation.frequencies.nu(nextpointidx, next_freq_idx) * next_shift
                - model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx) * curr_shift;
            const Real dtau = delta_tau(rayposidx, next_freq_idx);

            const Real impl_term = (dtau + expm1(-dtau)) / dtau;
            const Real impl_freq_der =
                dIdnu_coef2_next(rayposidx, next_freq_idx)
                    * intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx))
                + dIdnu_coef3_next(rayposidx, next_freq_idx)
                      * intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx));

            // Do the implicit part of the computation: add current intensity and source term, then
            // divide by the (1-freq derivative) term
            intensities(rayposidx, next_freq_idx) =
                (intensities(rayposidx, next_freq_idx)
                    + impl_term * S_next(rayposidx, next_freq_idx) // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - dIdnu_coef1_next(rayposidx, next_freq_idx) * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    }
}

///  Solver for the non-approximate comoving shortcharacteristics method. This is a first order in
///  space, second order in frequency accurate solver.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] o: point index of the origin of the ray
///   @param[in] r: ray direction index
///   @param[in] rayidx: index of ray traced through the model in the direction
///   @param[in] dshift_max: maximum doppler shift (currently does not do anything)
///   @note This solver might suffer a bit from spatial inaccuracies, as the traced rays might not
///   go exactly through the positions
template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_comoving_order_2_sparse(Model& model,
    const Size o,      // ray origin point
    const Size r,      // ray direction index
    const Size rayidx, // ray trace index
    const double dshift_max) {
    const Size rr = r;
    const Size ar = model.geometry.rays.antipod[rr];

    Matrix<Real>& intensities = intensities_(); // will contain all computed intensities

    Vector<Size>& nr = nr_();

    Vector<double>& shift =
        shift_(); // Contains the doppler shifts for all points on the ray. This
                  // need to be changed to 2.0-shift to get the shift direction
                  // I want for my frequencies for the forward ray..
                  // This is due to the default computed shift mapping from
                  // static frame to comoving frame (and I want to do the reverse).

    double Z  = 0.0; // distance along ray
    double dZ = 0.0; // last distance increment

    Size first_interesting_rayposidx = centre;
    Size last_interesting_rayposidx  = centre;

    // Trace the ray, getting all points on the ray and indicating until which point I actually need
    // to compute anything. Note: technically, it does not matter which frame I use, as long as it
    // is a fixed frame for the entire ray
    first_() = trace_ray_comoving<Rest, use_adaptive_directions>(model.geometry, o, rr, rr, rayidx, dshift_max, -1,
                   centre - 1, centre - 1, first_interesting_rayposidx)
             + 1;
    last_() = trace_ray_comoving<Rest, use_adaptive_directions>(model.geometry, o, ar, rr, rayidx, dshift_max, +1,
                  centre + 1, centre, last_interesting_rayposidx)
            - 1;

    nr_()[centre]    = o;
    shift_()[centre] = model.geometry.get_shift<Rest, use_adaptive_directions>(o, rr, o, 0);
    n_tot_()         = (last_() + 1) - first_();

    // Now is the perfect time to setup the boundary conditions and data for the forward ray
    comoving_ray_bdy_setup_forward<approx>(model, last_interesting_rayposidx);

    if (intensity_origin[nr[first_()]].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nr[first_()]);
        Size rcur  = std::get<1>(valid_rcur);
    // Check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the
    // if-clause
    // if (closest_ray(rr, nr[first_()]) == rayidx) {
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nr[first_()], rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nr[first_()], freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nr[first_()], k) +=
                lspec.quadrature.weights[z] * wt * intensities(first_(), freqid);
            // Lambda term does not apply to boundary points
        }
    }

    Size rayposidx = first_() + 1; // ray position index -> point index through nr[rayposidx]
    // from first_+1 to last_ index, trace the ray in
    while (rayposidx <= last_interesting_rayposidx) {
        const Real shift_next     = 2.0 - shift_()[rayposidx];
        const Real shift_curr     = 2.0 - shift_()[rayposidx - 1];
        const bool is_upward_disc = (shift_next >= shift_curr);

        // solve_comoving_single_step(model, rayposidx, rayidx, r, is_upward_disc, true);
        solve_comoving_single_step<use_adaptive_directions>(model, rayposidx, o, r, is_upward_disc, true);
        rayposidx++;
    }

    // Now is the perfect time to setup the boundary conditions and data for the backward ray
    comoving_ray_bdy_setup_backward<approx>(model, first_interesting_rayposidx);

    if (intensity_origin[nr[last_()]].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nr[last_()]);
        Size rcur  = std::get<1>(valid_rcur);
    // Check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the
    // if-clause
    // if (closest_ray(rr, nr[last_()]) == rayidx) {
        // std::cout<<"setting bdy intensities at last"<<std::endl;
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nr[last_()], rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nr[last_()], freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nr[last_()], k) +=
                lspec.quadrature.weights[z] * wt * intensities(last_(), freqid);
            // Lambda term does not apply to boundary points
        }
    }

    rayposidx = last_() - 1; // ray position index -> point index through nr[rayposidx]
    // From last_()-1 to first_ index, we trace the ray backwards; Warning: rayposidx is an unsigned
    // int, so +1 necessary to both sides
    while (rayposidx + 1 >= first_interesting_rayposidx + 1) {
        const Real shift_next     = shift_()[rayposidx];
        const Real shift_curr     = shift_()[rayposidx + 1];
        const bool is_upward_disc = (shift_next >= shift_curr);

        // solve_comoving_single_step(model, rayposidx, rayidx, r, is_upward_disc, false);
        solve_comoving_single_step<use_adaptive_directions>(model, rayposidx, o, r, is_upward_disc, false);
        rayposidx--;
    }
}

///  Solver for the approximate comoving shortcharacteristics method. This is a first order in
///  space, second order in frequency accurate solver. Should only be used when no significant
///  non-monotonic velocity fields are present.
///   @param[in] model: the model to apply the comoving solver to
///   @note This solver might suffer a bit from spatial inaccuracies, as the traced rays might not
///   go exactly through the positions
template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_comoving_local_approx_order_2_sparse(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        accelerated_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k) = 0.0;
            }
        })
    }

    std::cout << "hnrays: " << model.parameters->hnrays() << std::endl;

    std::vector<Size> cum_n_points_to_trace_ray_through(model.parameters->hnrays() + 1, 0);
    Size counter = 0; // will contain the total number of points to trace through
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
        counter += points_to_trace_ray_through[rr].size();
        cum_n_points_to_trace_ray_through[rr + 1] = counter;
    }

    // Parallelization over both rays and ray directions, as otherwise the load balancing might be a
    // bit iffy for smaller models
    accelerated_for(totalidx, counter, {
        auto cum_n_points_pointer = std::upper_bound(cum_n_points_to_trace_ray_through.begin(),
            cum_n_points_to_trace_ray_through.end(), totalidx);
        Size rr =
            std::distance(cum_n_points_to_trace_ray_through.begin(), cum_n_points_pointer) - 1;
        Size rayidx             = totalidx - cum_n_points_to_trace_ray_through[rr];
        const Size o            = points_to_trace_ray_through[rr][rayidx];
        const double dshift_max = get_dshift_max(model, o);
        solve_comoving_local_approx_order_2_sparse<approx, use_adaptive_directions>(model, o, rr, rayidx, dshift_max);
    })
}

///  Solver for the approximate comoving shortcharacteristics method. This is a first order in
///  space, second order in frequency accurate solver. Should only be used when no significant
///  non-monotonic velocity fields are present.
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] o: point index of the origin of the ray
///   @param[in] r: ray direction index
///   @param[in] rayidx: index of ray traced through the model in the direction
///   @param[in] dshift_max: maximum doppler shift (currently does not do anything)
///   @note This solver might suffer a bit from spatial inaccuracies, as the traced rays might not
///   go exactly through the positions
template <ApproximationType approx, bool use_adaptive_directions>
accel inline void Solver ::solve_comoving_local_approx_order_2_sparse(Model& model,
    const Size o,      // ray origin point
    const Size r,      // ray direction index
    const Size rayidx, // ray trace index
    const double dshift_max) {

    const Size rr = r;
    const Size ar = model.geometry.rays.antipod[rr];

    Vector<Real>& cma_intensities =
        cma_computed_intensities_(); // will computed intensities for a single step

    Vector<Size>& nr = nr_();

    Vector<double>& shift =
        shift_(); // Contains the doppler shifts for all points on the ray. For the forward ray, we
                  // need to do 2.0-shift to get the shift direction I want for my frequencies.
                  // This is due to the default computed shift mapping from static frame to comoving
                  // frame (I want to do the reverse) for the backward ray, this should be correct
                  // (as we need to do the reverse)

    double Z  = 0.0; // distance along ray
    double dZ = 0.0; // last distance increment

    // Trace the ray, getting all points on the ray and checking until which point we actually need
    // to use the comoving solver Note: technically, it does not matter which frame I use, as long
    // as it is a fixed frame for the entire ray
    Size first_interesting_rayposidx = centre;
    Size last_interesting_rayposidx  = centre;
    first_() = trace_ray_comoving<Rest, use_adaptive_directions>(model.geometry, o, rr, rr, rayidx, dshift_max, -1,
                   centre - 1, centre - 1, first_interesting_rayposidx)
             + 1;
    last_() = trace_ray_comoving<Rest, use_adaptive_directions>(model.geometry, o, ar, rr, rayidx, dshift_max, +1,
                  centre + 1, centre, last_interesting_rayposidx)
            - 1;

    nr_()[centre]    = o;
    shift_()[centre] = model.geometry.get_shift<Rest, use_adaptive_directions>(o, rr, o, 0);
    n_tot_()         = (last_() + 1) - first_();

    // Doppler shifts are computed in the opposite direction as usual, so correcting for this
    for (Size temp_rayidx = first_(); temp_rayidx <= last_(); temp_rayidx++) {
        shift_()[temp_rayidx] = 2.0 - shift_()[temp_rayidx];
    }

    // and also initialize helper variable for solver
    Vector<char>& cma_compute_curr_opacity = cma_compute_curr_opacity_();
    // also compute boundary intensities to start
    for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
        cma_compute_curr_opacity[freqidx] = true;

        cma_intensities[freqidx] = boundary_intensity(model, nr_()[first_()],
            model.radiation.frequencies.sorted_nu(nr_()[first_()], freqidx));
    }

    if (intensity_origin[nr[first_()]].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nr[first_()]);
        Size rcur  = std::get<1>(valid_rcur);
    // check if closest ray of the starting boundary point // maybe todo: replace with some weights
    // 0/1 for eliminating the if-clause
    // if (closest_ray(rr, nr[first_()]) == rayidx) {
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nr[first_()], rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nr[first_()], freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nr[first_()], k) += lspec.quadrature.weights[z] * wt * cma_intensities[freqid];
            // Lambda term does not apply to the boundary
        }
    }

    Size rayposidx = first_() + 1; // ray position index -> point index through nr[rayposidx]
    // from first_+1 to the last interesting index, trace the ray forwards
    while (rayposidx <= last_interesting_rayposidx) {
        const Real shift_next     = shift_()[rayposidx];
        const Real shift_curr     = shift_()[rayposidx - 1];
        const bool is_upward_disc = (shift_next >= shift_curr);

        // Mapping and solving can be done for position increment individually, because the
        // approximate boundary conditions do not retain memory of intensities at previous positions
        comoving_local_approx_map_data<approx>(model, nr_()[rayposidx - 1], nr_()[rayposidx],
            shift_curr, shift_next, is_upward_disc, dZ_()[rayposidx - 1], nr_()[first_()]);
        // solve_comoving_local_approx_single_step(model, nr_()[rayposidx], rayidx, r, is_upward_disc);
        solve_comoving_local_approx_single_step<use_adaptive_directions>(model, nr_()[rayposidx], o, r, is_upward_disc);
        rayposidx++;
    }

    // Doppler shifts in opposite direction evidently need to have the opposite shift
    for (Size temp_rayidx = first_(); temp_rayidx <= last_(); temp_rayidx++) {
        shift_()[temp_rayidx] = 2.0 - shift_()[temp_rayidx];
    }

    // Compute boundary intensities and setup helper variable
    for (Size freqidx = 0; freqidx < model.parameters->nfreqs(); freqidx++) {
        cma_compute_curr_opacity[freqidx] = true;

        cma_intensities[freqidx] = boundary_intensity(
            model, nr_()[last_()], model.radiation.frequencies.sorted_nu(nr_()[last_()], freqidx));
    }

    if (intensity_origin[nr[last_()]].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nr[last_()]);
        Size rcur  = std::get<1>(valid_rcur);
    // check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the
    // if-clause
    // if (closest_ray(rr, nr[last_()]) == rayidx) {
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nr[last_()], rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nr[last_()], freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nr[last_()], k) += lspec.quadrature.weights[z] * wt * cma_intensities[freqid];
            // Lambda term does not apply to the boundary
        }
    }

    rayposidx = last_() - 1; // ray position index -> point index through nr[rayposidx]
    // From last_()-1 to the first interesting index, we trace the ray backwards; Warning: rayposidx
    // is an unsigned int, so +1 necessary to both sides
    while (rayposidx + 1 >= first_interesting_rayposidx + 1) {
        const Real shift_next     = shift_()[rayposidx];
        const Real shift_curr     = shift_()[rayposidx + 1];
        const bool is_upward_disc = (shift_next >= shift_curr);

        // Mapping and solving can be done for position increment individually, because the
        // approximate boundary conditions do not retain memory of intensities at previous positions
        comoving_local_approx_map_data<approx>(model, nr_()[rayposidx + 1], nr_()[rayposidx],
            shift_curr, shift_next, is_upward_disc, dZ_()[rayposidx], nr_()[last_()]);
        // solve_comoving_local_approx_single_step(model, nr_()[rayposidx], rayidx, r, is_upward_disc);
        solve_comoving_local_approx_single_step<use_adaptive_directions>(model, nr_()[rayposidx], o, r, is_upward_disc);
        rayposidx--;
    }
}

///  Approximate comoving solver helper function: Matches the sorted frequency indices such that the
///  frequency difference is minimal.
///   @param[in] model: to model to apply the approximate comoving solver to
///   @param[in] curr_point: point index of the current point
///   @param[in] next_point: point index of the next point
///   @param[in] curr_shift: doppler shift of the current point
///   @param[in] next_shift: doppler shift of the next point
///   @param[in] is_upward_disc: determines where to place the frequency boundary conditions
///   @param[in] dZ: position increment
///   @param[in] bdy_point: point index of the start of the ray
template <ApproximationType approx>
accel inline void Solver ::comoving_local_approx_map_data(Model& model, const Size curr_point,
    const Size next_point, const Real curr_shift, const Real next_shift, const bool is_upward_disc,
    const Real dZ, const Size bdy_point) {
    bool is_in_bounds = false;

    Size in_bounds_count =
        0; // to decide whether to artificially add some boundary to the implicit part (if < 2) then
           // the implicit freq derivative cannot be computed up to second order
    Size curr_freq_count = 0; // count for explicit part to determine when to say that something is
                              // outside of our bounds
    std::set<Size> encountered_lines = std::set<Size>();

    if (is_upward_disc) {
        // Starting from the highest frequency
        //+1 to all indices due to using unsigned ints when looping down (overflow otherwise)
        Size curr_freq_idx = model.parameters->nfreqs() - 1 + 1;
        Size next_freq_idx = model.parameters->nfreqs() - 1 + 1;

        Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(
            curr_point, curr_freq_idx - 1); // arbitrary frame, as within a single point, one does
                                            // not need to care about doppler shifts
        Size curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        encountered_lines.insert(curr_line_idx);

        const Size min_freq_idx = 0 + 1;

        // Match frequency indices
        while ((curr_freq_idx >= (min_freq_idx + 1)) && next_freq_idx >= min_freq_idx) {
            // check if the frequency at the previous point index is higher than this frequency (in
            // static frame)
            if (model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx - 1) * curr_shift
                > model.radiation.frequencies.sorted_nu(next_point, next_freq_idx - 1)
                      * next_shift) {
                // the frequency of the current point at curr_freq_idx is too high, thus the index
                // needs to be lowered.
                curr_freq_idx--;

                unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(curr_point,
                    curr_freq_idx - 1); // arbitrary frame, as within a single point, one does not
                                        // need to care about doppler shifts
                curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
                encountered_lines.insert(curr_line_idx);

                // whenever we encounter a new line, the size of encountered_lines will increase
                // This condition is only satisfied whenever we encounter a new line after all quads
                // of all previous lines are encountered. Thus we define this to be true when
                // encountering a gap in the current (explicit) frequency quadrature
                if (curr_freq_count
                    == (encountered_lines.size() - 1) * model.parameters->nquads()) {
                    is_in_bounds    = false;
                    in_bounds_count = 0;
                } else {
                    is_in_bounds = true;
                }
                curr_freq_count++; // increment after checking, otherwise off by one error

            } else // freq matching index must be determined if it is just higher than the freq at
                   // the previous point
            {
                const Real curr_freq =
                    model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx - 1)
                    * curr_shift; // index - 1 due to size
                const Real next_freq =
                    model.radiation.frequencies.sorted_nu(next_point, next_freq_idx - 1)
                    * next_shift; // index - 1 due to size

                is_in_bounds = (curr_freq_count
                                > (encountered_lines.size() - 1) * model.parameters->nquads() + 2);
                comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ,
                    is_in_bounds, in_bounds_count, next_freq_idx - 1, curr_freq_idx - 1, next_freq,
                    curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

                next_freq_idx--;
                in_bounds_count++;
            }
        }

        // match remaining points as best as possible (with the outermost point)
        while (next_freq_idx >= min_freq_idx) {
            const Real curr_freq =
                model.radiation.frequencies.sorted_nu(curr_point, min_freq_idx - 1)
                * curr_shift; // index - 1 due to size
            const Real next_freq =
                model.radiation.frequencies.sorted_nu(next_point, next_freq_idx - 1)
                * next_shift; // index - 1 due to size
            is_in_bounds = (next_freq >= curr_freq);
            comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds,
                in_bounds_count, next_freq_idx - 1, min_freq_idx - 1, next_freq, curr_freq,
                next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

            next_freq_idx--;
            in_bounds_count++;
        }

    } else {
        // Starting from the lowest frequency
        Size curr_freq_idx      = 0;
        Size next_freq_idx      = 0;
        const Size max_freq_idx = model.parameters->nfreqs() - 1;

        Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(
            curr_point, curr_freq_idx); // arbitrary frame, as within a single point, one does not
                                        // need to care about doppler shifts
        Size curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        encountered_lines.insert(curr_line_idx);

        // Match frequency indices
        while ((curr_freq_idx <= (max_freq_idx - 1)) && next_freq_idx <= max_freq_idx) {
            // check if a previous frequency index exists which is higher than this (in static
            // frame)
            if (model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx) * curr_shift
                < model.radiation.frequencies.sorted_nu(next_point, next_freq_idx) * next_shift) {
                // the frequency of the current point at curr_freq_idx is too low, thus the index
                // needs to be incremented.
                curr_freq_idx++;

                unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(
                    curr_point, curr_freq_idx); // arbitrary frame, as within a single point, one
                                                // does not need to care about doppler shifts
                curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
                encountered_lines.insert(curr_line_idx);

                // whenever we encounter a new line, the size of encountered_lines will increase
                // This condition is only satisfied whenever we encounter a new line after all quads
                // of all previous lines are encountered Thus we define this to be true when
                // encountering a gap in the current (explicit) frequency quadrature
                if (curr_freq_count
                    == (encountered_lines.size() - 1) * model.parameters->nquads()) {
                    is_in_bounds    = false;
                    in_bounds_count = 0;
                } else {
                    is_in_bounds = true;
                }
                curr_freq_count++; // increment after checking, otherwise off by one error

            } else // freq matching index must be determined if it is just higher than the freq at
                   // the previous point
            {
                const Real curr_freq =
                    model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx) * curr_shift;
                const Real next_freq =
                    model.radiation.frequencies.sorted_nu(next_point, next_freq_idx) * next_shift;
                is_in_bounds = (curr_freq_count
                                > (encountered_lines.size() - 1) * model.parameters->nquads() + 2);
                comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ,
                    is_in_bounds, in_bounds_count, next_freq_idx, curr_freq_idx, next_freq,
                    curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

                next_freq_idx++;
                in_bounds_count++;
            }
        }

        // match remaining points as best as possible (with the outermost point)
        while (next_freq_idx <= max_freq_idx) {
            const Real curr_freq =
                model.radiation.frequencies.sorted_nu(curr_point, max_freq_idx) * curr_shift;
            const Real next_freq =
                model.radiation.frequencies.sorted_nu(next_point, next_freq_idx) * next_shift;
            is_in_bounds = (next_freq <= curr_freq);
            comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds,
                in_bounds_count, next_freq_idx, max_freq_idx, next_freq, curr_freq, next_shift,
                curr_shift, curr_line_idx, is_upward_disc, bdy_point);

            next_freq_idx++;
            in_bounds_count++;
        }
    }

    // also map the explicit frequency derivative now
    if (is_upward_disc) // determines the direction of the derivative to choose
    {
        // Iteration over next_freq_idx, as we want to compute intensities at the next point at
        // those frequencies. The corresponding frequencies at the current point will just be
        // derived from the frequency matching.
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];
            if (curr_freq_idx
                >= model.parameters->nfreqs()
                       - 2) // making sure that we can compute the explicit second order forward
                            // frequency derivative (as we need two extra points in that direction)
            {
                cma_dIdnu_expl_()[next_freq_idx] =
                    0.0;  // not enough data, so assume the derivative to be zero (in a slightly
                          // better implementation, we could use a first order derivative for the
                          // second to last frequency index)
                continue; // next iteration might still contain some interesting stuff to compute
            }
            const Real curr_freq =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx) * curr_shift;
            const Real curr_freqp1 =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx + 1) * curr_shift;
            const Real curr_freqp2 =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx + 2) * curr_shift;
            const Real curr_I   = cma_computed_intensities_()[curr_freq_idx];
            const Real curr_Ip1 = cma_computed_intensities_()[curr_freq_idx + 1];
            const Real curr_Ip2 = cma_computed_intensities_()[curr_freq_idx + 2];

            const Real dfreqsmall =
                curr_freqp2 - curr_freqp1; // cma_start_frequencies_()[next_freq_idx+2] -
                                           // cma_start_frequencies_()[next_freq_idx+1];
            const Real dfreqlarge =
                curr_freqp2 - curr_freq; // cma_start_frequencies_()[next_freq_idx+2] -
                                         // cma_start_frequencies_()[next_freq_idx+0];

            const Real dIdnu_coef3_curr =
                -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
            const Real dIdnu_coef2_curr =
                dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
            const Real dIdnu_coef1_curr = -dIdnu_coef3_curr - dIdnu_coef2_curr; // curr point itself

            const Real dIdnu_expl = curr_I * dIdnu_coef1_curr + curr_Ip1 * dIdnu_coef2_curr
                                  + curr_Ip2 * dIdnu_coef3_curr;
            cma_dIdnu_expl_()[next_freq_idx] = dIdnu_expl;
        }
    } else {
        // Iteration over next_freq_idx, as we want to compute intensities at the next point at
        // those frequencies. The corresponding frequencies at the current point will just be
        // derived from the frequency matching.
        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];
            if (curr_freq_idx
                <= 1) // making sure that we can compute the explicit second order forward frequency
                      // derivative (as we need two extra points in that direction)
            {
                cma_dIdnu_expl_()[next_freq_idx] =
                    0.0;  // not enough data, so assume the derivative to be zero (in a slightly
                          // better implementation, we could use a first order derivative for the
                          // second frequency index)
                continue; // next iteration might still contain some interesting stuff to compute
            }
            const Real curr_freq =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx) * curr_shift;
            const Real curr_freqm1 =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx - 1) * curr_shift;
            const Real curr_freqm2 =
                model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx - 2) * curr_shift;
            const Real curr_I   = cma_computed_intensities_()[curr_freq_idx];
            const Real curr_Im1 = cma_computed_intensities_()[curr_freq_idx - 1];
            const Real curr_Im2 = cma_computed_intensities_()[curr_freq_idx - 2];

            const Real dfreqsmall =
                curr_freqm2 - curr_freqm1; // cma_start_frequencies_()[next_freq_idx+2] -
                                           // cma_start_frequencies_()[next_freq_idx+1];
            const Real dfreqlarge =
                curr_freqm2 - curr_freq; // cma_start_frequencies_()[next_freq_idx+2] -
                                         // cma_start_frequencies_()[next_freq_idx+0];

            const Real dIdnu_coef3_curr =
                -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
            const Real dIdnu_coef2_curr =
                dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
            const Real dIdnu_coef1_curr = -dIdnu_coef3_curr - dIdnu_coef2_curr; // curr point itself

            const Real dIdnu_expl = curr_I * dIdnu_coef1_curr + curr_Im1 * dIdnu_coef2_curr
                                  + curr_Im2 * dIdnu_coef3_curr;
            cma_dIdnu_expl_()[next_freq_idx] = dIdnu_expl;
        }
    }
}

///  Approximate comoving solver helper function: map a single datapoint for the approximate
///  comoving solver
///   @param[in] model: the model to apply the approximate comoving solver to
///   @param[in] curr_point: the point index of the current point
///   @param[in] next_point: the point index of the next point
///   @param[in] dZ: the distance increment
///   @param[in] is_in_bounds: whether or not the next frequency should be computed using
///   (approximate) boundary conditions
///   @param[in] curr_freq_count: is the nth frequency we try to map within a consecutive line
///   region. Determines whether to use (approximate) boundary conditions
///   @param[in] next_freq_idx: frequency index of the mapped frequency at the next point
///   @param[in] curr_freq_idx: frequency index of the mapped frequency at the current point
///   @param[in] next_freq: mapped frequency in comoving frame at next point
///   @param[in] curr_freq: mapped frequency in comoving frame at current point
///   @param[in] next_shift: doppler shift at next point
///   @param[in] curr_shift: doppler shift at current point
///   @param[in] curr_line_idx: line index of the current line (only for OneLine approx)
///   @param[in] is_upward_disc: determines the ordering of the implicit part of the computations
///   @param[in] bdy_point: point index of the boundary point. Used for determining the boundary
///   intensities
///   @warning Assumes the frequencies to be set in a fixed order, depending on is_upward_disc. See
///   code itself for more details.
template <ApproximationType approx>
accel inline void Solver ::comoving_approx_map_single_data(Model& model, const Size curr_point,
    const Size next_point, const Real dZ, const bool is_in_bounds, const Size curr_freq_count,
    const Size next_freq_idx, const Size curr_freq_idx, const Real next_freq, const Real curr_freq,
    const Real next_shift, const Real curr_shift, const Size curr_line_idx,
    const bool is_upward_disc, const Size bdy_point) {
    // Note: data mapping is done using next_freq_idx (as we want to compute intensities at all
    // those frequencies)

    // If boundary conditions are required, then use the planck intensity
    if (!is_in_bounds || curr_freq_count < 2) {

        // Define new start intensities using boundary conditions
        cma_start_intensities_()[next_freq_idx] =
            boundary_intensity(model, bdy_point, next_freq / curr_shift);

        cma_start_frequencies_()[next_freq_idx] = next_freq;
        cma_start_frequency_index_()[next_freq_idx] =
            is_upward_disc ? (model.parameters->nfreqs() - 1)
                           : 0; // dummy behavior, signaling that this is a bdy point for the
                                // different directions, a different frequency is chosen for
                                // simplicity (otherwise I need another 'is_bdy' field)
        cma_end_frequencies_()[next_freq_idx] = next_freq;

        const Real shift_c = curr_shift;
        const Real shift_n = next_shift;
        // make sure that the computed intensity will not change significantly from the start
        // intensity
        Real Snext                = 0.0;
        Real Scurr                = 0.0;
        Real dtau                 = model.parameters->comoving_min_dtau;
        Real chi_next             = model.parameters->min_opacity;
        Real chi_curr             = 0.0;  // As we are out of bounds either way
        bool compute_curr_opacity = true; // we have no information at the bdy freq index, so we
                                          // formally need to compute the opacity; mainly important
                                          // if we ever get around to implementing continuum stuff

        // as usual, I define the shift backwards to obtain observer frame frequencies once again
        compute_source_dtau<approx>(model, curr_point, next_point, curr_line_idx,
            next_freq / shift_c, next_freq / shift_n, shift_c, shift_n, dZ, compute_curr_opacity,
            dtau, chi_curr, chi_next, Scurr, Snext);
        // cheat on dIdnu coeffs, as these are all zero on the boundary
        cma_dIdnu_coef1_next_()[next_freq_idx] = 0;
        cma_dIdnu_coef2_next_()[next_freq_idx] = 0;
        cma_dIdnu_coef3_next_()[next_freq_idx] = 0;

        cma_chi_next_()[next_freq_idx]             = chi_next;
        cma_compute_next_opacity_()[next_freq_idx] = compute_curr_opacity; // for the next thing

        cma_S_next_()[next_freq_idx] = Snext;
        cma_S_curr_()[next_freq_idx] = Scurr;
        // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
        dtau                            = std::max(dtau, model.parameters->comoving_min_dtau);
        cma_delta_tau_()[next_freq_idx] = dtau;
    } else {
        // Set starting intensities, frequencies to the correct computed intensity of the previous
        // iteration
        cma_start_intensities_()[next_freq_idx]     = cma_computed_intensities_()[curr_freq_idx];
        cma_start_frequencies_()[next_freq_idx]     = curr_freq;
        cma_start_frequency_index_()[next_freq_idx] = curr_freq_idx;
        cma_end_frequencies_()[next_freq_idx]       = next_freq;

        const Real shift_c = curr_shift;
        const Real shift_n = next_shift;

        // Set some dummy values which will be overwritten
        Real Snext    = 0.0;
        Real dtau     = model.parameters->comoving_min_dtau;
        Real chi_next = model.parameters->min_opacity;
        Real Scurr = cma_S_curr_()[curr_freq_idx]; // if we have already computed the opacity in the
                                                   // low doppler shift case, the current source
                                                   // will not be computed...

        Real chi_curr = cma_chi_curr_()[curr_freq_idx]; // retrieving the correct frequency
        bool compute_curr_opacity = cma_compute_curr_opacity_()[curr_freq_idx];

        compute_source_dtau<approx>(model, curr_point, next_point, curr_line_idx,
            curr_freq / shift_c, next_freq / shift_n, shift_c, shift_n, dZ, compute_curr_opacity,
            dtau, chi_curr, chi_next, Scurr, Snext);

        cma_S_next_()[next_freq_idx] = Snext;
        cma_S_curr_()[next_freq_idx] = Scurr;

        cma_chi_next_()[next_freq_idx] = chi_next;
        cma_chi_curr_()[curr_freq_idx] = chi_curr;
        cma_compute_curr_opacity_()[curr_freq_idx] =
            compute_curr_opacity; // if we somehow encounter a current frequency index twice, we set
                                  // this
        cma_compute_next_opacity_()[next_freq_idx] = compute_curr_opacity; // for the next thing

        // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
        dtau                            = std::max(dtau, model.parameters->comoving_min_dtau);
        cma_delta_tau_()[next_freq_idx] = dtau;

        // WARNING: Assumes the frequencies to be set in a fixed order!
        // if is_upward_disc==TRUE, then the frequency indices should be traversed from nfreqs-1 to
        // 0 if is_upward_disc==FALSE, then the frequency indices should be traversed from 0 to
        // nfreqs-1. This is due to cma_end_frequencies_ only being set one at a time.
        if (is_upward_disc) {
            // And now do the same (but simpler; less index management) for the implicit part
            Real dfreqsmall = cma_end_frequencies_()[next_freq_idx + 2]
                            - cma_end_frequencies_()[next_freq_idx + 1];

            Real dfreqlarge = cma_end_frequencies_()[next_freq_idx + 2]
                            - cma_end_frequencies_()[next_freq_idx + 0];

            cma_dIdnu_coef3_next_()[next_freq_idx] =
                -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
            cma_dIdnu_coef2_next_()[next_freq_idx] =
                dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
            cma_dIdnu_coef1_next_()[next_freq_idx] =
                -cma_dIdnu_coef3_next_()[next_freq_idx]
                - cma_dIdnu_coef2_next_()[next_freq_idx]; // curr point itself;

        } else {
            // And now do the same (but simpler; less index management) for the implicit part
            // Note: minus signs in front, as we are now computing the first derivative using points
            // on the other side
            Real dfreqsmall = -(cma_end_frequencies_()[next_freq_idx - 1]
                                - cma_end_frequencies_()[next_freq_idx - 2]);

            Real dfreqlarge = -(cma_end_frequencies_()[next_freq_idx - 0]
                                - cma_end_frequencies_()[next_freq_idx - 2]);

            cma_dIdnu_coef3_next_()[next_freq_idx] =
                -dfreqsmall / (std::pow(dfreqlarge, 2.0) - dfreqlarge * dfreqsmall); // farthest
            cma_dIdnu_coef2_next_()[next_freq_idx] =
                dfreqlarge / (-std::pow(dfreqsmall, 2.0) + dfreqlarge * dfreqsmall); // nearer
            cma_dIdnu_coef1_next_()[next_freq_idx] =
                -cma_dIdnu_coef3_next_()[next_freq_idx]
                - cma_dIdnu_coef2_next_()[next_freq_idx]; // curr point itself;
        }
        // end of non-boundary data mapping
    }
}

///  Approximate comoving solver: solves the comoving equations for a single step
///   @param[in] model: the model to apply the comoving solver to
///   @param[in] next_point: point index of the next point
// DEPRECATED///   @param[in] rayidx: index of ray to trace in the given direction
///   @param[in] o: position index of ray origin used to trace ray
///   @param[in] rr: ray direction index to check whether to save the current intensities (if the
///   ray lies closest in the current direction to the point)
///   @param[in] is_upward_disc: determines on which side of the frequency spectrum to put extra
///   boundary conditions
///   @note: assume the data to be setup correctly using Solver::comoving_local_approx_map_data
template <bool use_adaptive_directions>
accel inline void Solver ::solve_comoving_local_approx_single_step(Model& model,
    const Size next_point, const Size o, const Size rr, const bool is_upward_disc) {
    Vector<Size>& nr      = nr_();    // stores the exact point indices
    Vector<double>& shift = shift_(); // stores the shift; should already be in the correct
                                      // direction
    Vector<Real>& cma_delta_tau = cma_delta_tau_(); // stores the optical depths
    Vector<Real>& cma_S_curr    = cma_S_curr_();    // source function at current point
    Vector<Real>& cma_S_next    = cma_S_next_();    // source function at next point
    Vector<Real>& cma_computed_intensities =
        cma_computed_intensities_(); // storage for computed intensities (at next point)
    Vector<Real>& cma_start_intensities =
        cma_start_intensities_(); // contains the start intensities for the computatino
    Vector<Real>& cma_start_frequencies =
        cma_start_frequencies_(); // contains the frequencies at the current point
    Vector<Real>& cma_end_frequencies =
        cma_end_frequencies_(); // contains the frequencies at the next point

    Vector<Real>& cma_dIdnu_expl = cma_dIdnu_expl_(); // contains the computed explicit intensity
                                                      // derivatives with respect to the frequency

    // Frequency derivative coefficients
    Vector<Real>& cma_dIdnu_coef1_next = cma_dIdnu_coef1_next_();
    Vector<Real>& cma_dIdnu_coef2_next = cma_dIdnu_coef2_next_();
    Vector<Real>& cma_dIdnu_coef3_next = cma_dIdnu_coef3_next_();

    const Size nextpointidx = next_point;

    // Do the explicit part of the computation
    // discretization determines which points are most definitely boundary points (for which the
    // computation is simpler)
    if (is_upward_disc) {
        // There exist two boundary conditions at model.parameters->nfreqs()-1,-2

        // outermost boundary
        const Real bdy_dtau1      = cma_delta_tau[model.parameters->nfreqs() - 1];
        const Real bdy_expl_term1 = (-expm1(-bdy_dtau1) - bdy_dtau1 * exp(-bdy_dtau1)) / bdy_dtau1;
        cma_computed_intensities[model.parameters->nfreqs() - 1] =
            cma_start_intensities[model.parameters->nfreqs() - 1] * exp(-bdy_dtau1)
            + bdy_expl_term1
                  * cma_S_curr[model.parameters->nfreqs()
                               - 1]; // source term //no frequency derivative term, as the frequency
                                     // shift of the boundary is by definition 0
        // second boundary
        const Real bdy_dtau2      = cma_delta_tau[model.parameters->nfreqs() - 2];
        const Real bdy_expl_term2 = (-expm1(-bdy_dtau2) - bdy_dtau2 * exp(-bdy_dtau2)) / bdy_dtau2;
        cma_computed_intensities[model.parameters->nfreqs() - 2] =
            cma_start_intensities[model.parameters->nfreqs() - 2] * exp(-bdy_dtau2)
            + bdy_expl_term2
                  * cma_S_curr[model.parameters->nfreqs()
                               - 2]; // source term //no frequency derivative term, as the frequency
                                     // shift of the boundary is by definition 0

        for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs() - 2;
             next_freq_idx++) {

            const Real deltanu =
                cma_end_frequencies[next_freq_idx] - cma_start_frequencies[next_freq_idx];
            const Real dtau = cma_delta_tau[next_freq_idx];

            const Real expl_term     = (-expm1(-dtau) - dtau * exp(-dtau)) / dtau;
            const Real expl_freq_der = cma_dIdnu_expl[next_freq_idx];

            cma_computed_intensities[next_freq_idx] =
                cma_start_intensities[next_freq_idx] * exp(-dtau)
                + expl_term * cma_S_curr[next_freq_idx]       // source term
                + expl_term * expl_freq_der * deltanu / dtau; // frequency derivative term
        }
    } else { //! is_upward_disc
        // There exist two boundary conditions at 0, 1

        // outermost boundary
        const Real bdy_dtau1      = cma_delta_tau[0];
        const Real bdy_expl_term1 = (-expm1(-bdy_dtau1) - bdy_dtau1 * exp(-bdy_dtau1)) / bdy_dtau1;
        cma_computed_intensities[0] =
            cma_start_intensities[0] * exp(-bdy_dtau1)
            + bdy_expl_term1 * cma_S_curr[0]; // source term //no frequency derivative term, as the
                                              // frequency shift of the boundary is by definition 0
        // second boundary
        const Real bdy_dtau2      = cma_delta_tau[1];
        const Real bdy_expl_term2 = (-expm1(-bdy_dtau2) - bdy_dtau2 * exp(-bdy_dtau2)) / bdy_dtau2;
        cma_computed_intensities[1] =
            cma_start_intensities[1] * exp(-bdy_dtau2)
            + bdy_expl_term2 * cma_S_curr[1]; // source term //no frequency derivative term, as the
                                              // frequency shift of the boundary is by definition 0

        for (Size next_freq_idx = 2; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {

            const Real deltanu =
                cma_end_frequencies[next_freq_idx] - cma_start_frequencies[next_freq_idx];
            const Real dtau = cma_delta_tau[next_freq_idx];

            const Real expl_term     = (-expm1(-dtau) - dtau * exp(-dtau)) / dtau;
            const Real expl_freq_der = cma_dIdnu_expl[next_freq_idx];

            cma_computed_intensities[next_freq_idx] =
                cma_start_intensities[next_freq_idx] * exp(-dtau)
                + expl_term * cma_S_curr[next_freq_idx]       // source term
                + expl_term * expl_freq_der * deltanu / dtau; // frequency derivative term
        }
    }

    // Implicit part ordering depends on the discretization direction
    if (is_upward_disc) {
        // There exist two boundary conditions at model.parameters->nfreqs()-1,-2
        const Real bdy_dtau1  = cma_delta_tau[model.parameters->nfreqs() - 1];
        const Real impl_term1 = (bdy_dtau1 + expm1(-bdy_dtau1)) / bdy_dtau1;
        cma_computed_intensities[model.parameters->nfreqs() - 1] =
            cma_computed_intensities[model.parameters->nfreqs() - 1]
            + impl_term1 * cma_S_next[model.parameters->nfreqs() - 1]; // source term
        // no freq derivative term

        const Real bdy_dtau2  = cma_delta_tau[model.parameters->nfreqs() - 2];
        const Real impl_term2 = (bdy_dtau2 + expm1(-bdy_dtau2)) / bdy_dtau2;
        cma_computed_intensities[model.parameters->nfreqs() - 2] =
            cma_computed_intensities[model.parameters->nfreqs() - 2]
            + impl_term2 * cma_S_next[model.parameters->nfreqs() - 2]; // source term
        // no freq derivative term

        // In case of positive doppler shifts, we need to start from the largest frequencies (where
        // the extra boundary conditions lie). As we are using a reverse loop over unsigned
        // integers, the indices need to be incremented by 1.
        for (Size next_freq_idx = model.parameters->nfreqs() - 2; next_freq_idx > 0;
             next_freq_idx--) {

            const Real deltanu =
                cma_end_frequencies[next_freq_idx - 1] - cma_start_frequencies[next_freq_idx - 1];
            const Real dtau = cma_delta_tau[next_freq_idx - 1];

            const Real impl_term     = (dtau + expm1(-dtau)) / dtau;
            const Real impl_freq_der = cma_dIdnu_coef2_next[next_freq_idx - 1]
                                         * cma_computed_intensities[next_freq_idx + 1 - 1]
                                     + cma_dIdnu_coef3_next[next_freq_idx - 1]
                                           * cma_computed_intensities[next_freq_idx + 2 - 1];

            cma_computed_intensities[next_freq_idx - 1] =
                (cma_computed_intensities[next_freq_idx - 1]
                    + impl_term * cma_S_next[next_freq_idx - 1] // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - cma_dIdnu_coef1_next[next_freq_idx - 1] * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    } else {
        // There exist two boundary conditions at 0,1
        const Real bdy_dtau1  = cma_delta_tau[0];
        const Real impl_term1 = (bdy_dtau1 + expm1(-bdy_dtau1)) / bdy_dtau1;
        cma_computed_intensities[0] =
            cma_computed_intensities[0]
            + impl_term1 * cma_S_next[0]; // source term // no freq derivative term

        const Real bdy_dtau2  = cma_delta_tau[1];
        const Real impl_term2 = (bdy_dtau2 + expm1(-bdy_dtau2)) / bdy_dtau2;
        cma_computed_intensities[1] =
            cma_computed_intensities[1]
            + impl_term2 * cma_S_next[1]; // source term // no freq derivative term

        for (Size next_freq_idx = 2; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
            const Real deltanu =
                cma_end_frequencies[next_freq_idx] - cma_start_frequencies[next_freq_idx];
            const Real dtau = cma_delta_tau[next_freq_idx];

            const Real impl_term = (dtau + expm1(-dtau)) / dtau;
            const Real impl_freq_der =
                cma_dIdnu_coef2_next[next_freq_idx] * cma_computed_intensities[next_freq_idx - 1]
                + cma_dIdnu_coef3_next[next_freq_idx] * cma_computed_intensities[next_freq_idx - 2];

            cma_computed_intensities[next_freq_idx] =
                (cma_computed_intensities[next_freq_idx]
                    + impl_term * cma_S_next[next_freq_idx] // source term
                    + impl_term * impl_freq_der * deltanu / dtau)
                / (1.0
                    - cma_dIdnu_coef1_next[next_freq_idx] * impl_term * deltanu
                          / dtau); // freq derivative term
        }
    }

    if (intensity_origin[nextpointidx].count(std::tuple(o, rr)))
    {
        //We have assigned this ray to contain data for this nextpointindex, so I assume the ray direction to be valid
        std::tuple<bool, Size> valid_rcur = model.geometry.rays.get_correspoding_direction_index<use_adaptive_directions>(o, rr, nextpointidx);
        Size rcur  = std::get<1>(valid_rcur);
    // Finally increment J if this ray lies closest to the point in the raydirection rr
    // if (closest_ray(rr, nextpointidx) == rayidx) {
        const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(nextpointidx, rcur);
        // then obviously add (weighted) to J
        for (Size freqid = 0; freqid < model.parameters->nfreqs(); freqid++) {
            // Get the details about the line to which the current frequency belongs
            const Size unsorted_freqidx =
                model.radiation.frequencies.corresponding_nu_index(nextpointidx, freqid);
            const Size l = model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z = model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

            lspec.J(nextpointidx, k) += lspec.quadrature.weights[z] * wt
                                      * cma_computed_intensities[freqid]; // Su_()[centre];

            // Computing the ALI element
            const Real deltanu = cma_end_frequencies[freqid] - cma_start_frequencies[freqid];
            const Real dtau    = cma_delta_tau[freqid];
            const Real constant =
                lspec.linedata.A[k] * lspec.quadrature.weights[z]
                * model.geometry.rays
                      .weight[rr]; // for integrating over both frequency (weighted with profile
                                   // function) and angle. We also need the einstein A coefficient

            const Real source_term =
                (dtau + expm1(-dtau)); // '/dtau' has already been applied to simplify the fraction
            const Real lambdaterm = constant * source_term
                                  / (dtau - cma_dIdnu_coef1_next[freqid] * source_term * deltanu);

            lspec.lambda.add_element(nextpointidx, k, nextpointidx, lambdaterm);
        }
    }

    // Finally overwrite some variables with their value at the next point; in this way, we prepare
    // for the next step
    for (Size next_freq_idx = 0; next_freq_idx < model.parameters->nfreqs(); next_freq_idx++) {
        cma_S_curr[next_freq_idx]                  = cma_S_next[next_freq_idx];
        cma_chi_curr_()[next_freq_idx]             = cma_chi_next_()[next_freq_idx];
        cma_compute_curr_opacity_()[next_freq_idx] = cma_compute_next_opacity_()[next_freq_idx];
    }
}

// template <Frame frame>
inline Size Solver ::get_ray_lengths_max_new_imager(
    Model& model, Image& image, const Vector3D& ray_dir) {
    const Size npixels = image.ImX.size(); // is ImY.size() (number of pixels in image)

    std::vector<Size> ray_lengths;
    ray_lengths.resize(npixels); // TODO: later on, define which pixels need to
                                 // be used for raytracing (adaptive imager)/or
                                 // define new partial images (might be more
                                 // elegant)

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        ray_lengths[pixidx] =
            get_ray_length_new_imager(model.geometry, origin, start_bdy_point, ray_dir);
    });

    const Size max_ray_length = (*std::max_element(ray_lengths.begin(), ray_lengths.end()));
    return max_ray_length;
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_shortchar_order_0(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.I.resize(
            model.parameters->nrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize(model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();
    }

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar = model.geometry.rays.get_antipod_index(rr);
            // Approach which just accumulates the intensity
            // contributions as the ray is traced

            solve_shortchar_order_0<approx, use_adaptive_directions>(model, o, rr);
            solve_shortchar_order_0<approx, use_adaptive_directions>(model, o, ar);

            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                model.radiation.u(rr, o, f) =
                    0.5 * (model.radiation.I(rr, o, f) + model.radiation.I(ar, o, f));
                model.radiation.v(rr, o, f) =
                    0.5 * (model.radiation.I(rr, o, f) - model.radiation.I(ar, o, f));
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}

/// BUGGED: v computation is incorrect
template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_uv(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar = model.geometry.rays.get_antipod_index(rr);

            const Real dshift_max = get_dshift_max(model, o);

            nr_()[centre]    = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray<CoMoving, use_adaptive_directions>(
                           model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
                     + 1;
            last_() = trace_ray<CoMoving, use_adaptive_directions>(
                          model.geometry, o, ar, dshift_max, +1, centre + 1, centre)
                    - 1;
            n_tot_() = (last_() + 1) - first_();

            if (n_tot_() > 1) {
                for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                    solve_feautrier_order_2_uv<approx>(model, o, f);

                    model.radiation.u(rr, o, f) = Su_()[centre];
                    model.radiation.v(rr, o, f) = Sv_()[centre];
                }
            } else {
                for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                    model.radiation.u(rr, o, f) =
                        boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.v(rr, o, f) = 0.0;
                }
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.v.copy_ptr_to_vec();
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_sparse(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k) = 0.0;
            }
        })
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
            threaded_for(o, model.parameters->npoints(), {
                const Vector3D nn =
                    model.geometry.rays.get_direction<use_adaptive_directions>(o, rr);
                const Size ar = model.geometry.rays.get_antipod_index(rr);
                const Real wt =
                    model.geometry.rays.get_weight<use_adaptive_directions>(o, rr) * two;

                const Real dshift_max = get_dshift_max(model, o);

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray<CoMoving, use_adaptive_directions>(
                               model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
                         + 1;
                last_() = trace_ray<CoMoving, use_adaptive_directions>(
                              model.geometry, o, ar, dshift_max, +1, centre + 1, centre)
                        - 1;
                n_tot_() = (last_() + 1) - first_();

                if (n_tot_() > 1) {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            solve_feautrier_order_2<approx>(model, o, lspec.nr_line[o][k][z]);

                            lspec.J(o, k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda<approx, use_adaptive_directions>(
                                model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                } else {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            lspec.J(o, k) +=
                                lspec.quadrature.weights[z] * wt
                                * boundary_intensity(model, o,
                                    model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_anis(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_0.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k)       = 0.0;
                lspec.J2_0(o, k)    = 0.0;
                lspec.J2_1_Re(o, k) = 0.0;
                lspec.J2_1_Im(o, k) = 0.0;
                lspec.J2_2_Re(o, k) = 0.0;
                lspec.J2_2_Im(o, k) = 0.0;
            }
        })
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
            threaded_for(o, model.parameters->npoints(), {
                const Size ar = model.geometry.rays.get_antipod_index(rr);
                const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(o, rr);
                const Vector3D nn =
                    model.geometry.rays.get_direction<use_adaptive_directions>(o, rr);

                const Real wt_0    = inv_sqrt2 * (three * nn.z() * nn.z() - one);
                const Real wt_1_Re = -sqrt3 * nn.x() * nn.z();
                const Real wt_1_Im = -sqrt3 * nn.y() * nn.z();
                const Real wt_2_Re = half * sqrt3 * (nn.x() * nn.x() - nn.y() * nn.y());
                const Real wt_2_Im = sqrt3 * nn.x() * nn.y();

                const Real dshift_max = get_dshift_max(model, o);

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray<CoMoving, use_adaptive_directions>(
                               model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
                         + 1;
                last_() = trace_ray<CoMoving, use_adaptive_directions>(
                              model.geometry, o, ar, dshift_max, +1, centre + 1, centre)
                        - 1;
                n_tot_() = (last_() + 1) - first_();

                if (n_tot_() > 1) {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            solve_feautrier_order_2<approx>(model, o, lspec.nr_line[o][k][z]);

                            const Real du = lspec.quadrature.weights[z] * wt * Su_()[centre];

                            lspec.J(o, k) += two * du;
                            lspec.J2_0(o, k) += wt_0 * du;
                            lspec.J2_1_Re(o, k) += wt_1_Re * du;
                            lspec.J2_1_Im(o, k) += wt_1_Im * du;
                            lspec.J2_2_Re(o, k) += wt_2_Re * du;
                            lspec.J2_2_Im(o, k) += wt_2_Im * du;
                        }
                    }
                } else {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            const Real du =
                                lspec.quadrature.weights[z] * wt
                                * boundary_intensity(model, o,
                                    model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));

                            lspec.J(o, k) += two * du;
                            lspec.J2_0(o, k) += wt_0 * du;
                            lspec.J2_1_Re(o, k) += wt_1_Re * du;
                            lspec.J2_1_Im(o, k) += wt_1_Im * du;
                            lspec.J2_2_Re(o, k) += wt_2_Re * du;
                            lspec.J2_2_Im(o, k) += wt_2_Im * du;
                        }
                    }
                }
            })
        }
    }
}

// sparse Feautrier solver, but now does not put point without close lines on the ray
template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_sparse_pruned_rays(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k) = 0.0;
            }
        })
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
        const Size ar     = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
            threaded_for(o, model.parameters->npoints(), {
                const Real dshift_max = get_dshift_max(model, o);
                const Real wt     = model.geometry.rays.get_weight<use_adaptive_directions>(o, rr) * two;
                const Vector3D nn = model.geometry.rays.get_direction<use_adaptive_directions>(o, rr);

                // first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1,
                // centre-1) + 1; last_ () = trace_ray <CoMoving> (model.geometry, o, ar,
                // dshift_max, +1, centre+1, centre  ) - 1; n_tot_() = (last_()+1) - first_();

                // Basic ray tracing to check if ray length is at least 2 (so no direct boundary
                // condition required) By definition of the comoving frame, o defines the reference
                // frame. To correctly compute the intensity,
                //  we will at least include the origin of the ray and the neighbors in the pruning
                //  process.
                const Size tot_ray_length =
                    1 + model.geometry.get_ray_length<CoMoving, use_adaptive_directions>(o, rr, dshift_max)
                    + model.geometry.get_ray_length<CoMoving, use_adaptive_directions>(o, ar, dshift_max);

                // if (n_tot_() > 1)
                if (tot_ray_length > 1) {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        const Real line_frequency = lspec.linedata.frequency[k];
                        // Per line, trace which part we actually need to solve
                        // Is identical for all pruned rays
                        nr_()[centre]    = o;
                        shift_()[centre] = 1.0;

                        first_() = trace_ray_pruned<CoMoving, use_adaptive_directions>(model, o, rr, dshift_max, -1,
                                       centre - 1, centre - 1, line_frequency)
                                 + 1;
                        last_() = trace_ray_pruned<CoMoving, use_adaptive_directions>(model, o, ar, dshift_max, +1,
                                      centre + 1, centre, line_frequency)
                                - 1;
                        n_tot_() = (last_() + 1) - first_();
                        // std::cout<<"ntot: "<<n_tot_()<<std::endl;

                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            solve_feautrier_order_2<approx>(model, o, lspec.nr_line[o][k][z]);

                            lspec.J(o, k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda<approx, use_adaptive_directions>(model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                } else {
                    for (Size k = 0; k < lspec.linedata.nrad; k++) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            lspec.J(o, k) +=
                                lspec.quadrature.weights[z] * wt
                                * boundary_intensity(model, o,
                                    model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize(model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();
    }

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    distributed_for(rr, rr_loc, model.parameters->hnrays(),
        {
            cout << "--- rr = " << rr << endl;

            accelerated_for(o, model.parameters->npoints(), {
                const Size ar = model.geometry.rays.get_antipod_index(rr);

                const Real dshift_max = get_dshift_max(model, o);

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray<CoMoving, use_adaptive_directions>(
                               model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
                         + 1;
                last_() = trace_ray<CoMoving, use_adaptive_directions>(
                              model.geometry, o, ar, dshift_max, +1, centre + 1, centre)
                        - 1;
                n_tot_() = (last_() + 1) - first_();

                if (n_tot_() > 1) {
                    for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                        solve_feautrier_order_2<approx>(model, o, f);

                        model.radiation.u(rr_loc, o, f) = Su_()[centre];
                        model.radiation.J(o, f) +=
                            Su_()[centre] * two
                            * model.geometry.rays.get_weight<use_adaptive_directions>(o, rr);

                        update_Lambda<approx, use_adaptive_directions>(model, rr, f);
                    }
                } else {
                    for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                        model.radiation.u(rr_loc, o, f) =
                            boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                        model.radiation.J(o, f) +=
                            two * model.geometry.rays.get_weight<use_adaptive_directions>(o, rr)
                            * model.radiation.u(rr_loc, o, f);
                    }
                }
            })

            pc::accelerator::synchronize();
        })

        model.radiation.u.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();

    cout << "MPI gathering J..." << endl;
    model.radiation.MPI_reduce_J();
    cout << "Done MPI gathering J." << endl;

    cout << "MPI gathering Lambda..." << endl;
    // Gather contributions to the ALO
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.MPI_gather();
    }
    cout << "Done MPI gathering Lambda." << endl;
}


template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2(Model& model, const Size rr) {
    Image image = Image(model.geometry, model.radiation.frequencies, Intensity, rr);

    accelerated_for(o, model.parameters->npoints(), {
        const Size ar = model.geometry.rays.get_antipod_index(rr);

        const Real dshift_max = get_dshift_max(model, o);

        nr_()[centre]    = o;
        shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);
        ;

        first_() =
            trace_ray<Rest, false>(model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
            + 1;
        last_() =
            trace_ray<Rest, false>(model.geometry, o, ar, dshift_max, +1, centre + 1, centre) - 1;
        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_feautrier_order_2<approx>(model, o, f);

                image.I(o, f) = two * Su_()[last_()]
                              - boundary_intensity(
                                  model, nr_()[last_()], model.radiation.frequencies.nu(o, f));
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(o, f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2_new_imager(
    Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix) {
    Image image =
        Image(model.geometry, model.radiation.frequencies, Intensity, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future; (but then while loop will be
    // required anyway)
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);
        const Real dshift_max = get_dshift_max(model, closest_bdy_point);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);
        first_() = trace_ray_imaging<Rest>(model.geometry, origin, closest_bdy_point, ray_dir,
                       dshift_max, -1, Z, centre - 1, centre - 1)
                 + 1;
        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction

        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_feautrier_order_2<approx>(model, closest_bdy_point, f);
                image.I(pixidx, f) = two * Su_()[last_()]
                                   - boundary_intensity(model, nr_()[last_()],
                                       model.radiation.frequencies.nu(closest_bdy_point, f));
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = boundary_intensity(
                    model, closest_bdy_point, model.radiation.frequencies.nu(closest_bdy_point, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}


template <ApproximationType approx>
inline void Solver ::image_shortchar_order_0_new_imager(
    Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix) {
    Image image =
        Image(model.geometry, model.radiation.frequencies, Intensity, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future; (but then while loop will be
    // required anyway)
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);
        const Real dshift_max = get_dshift_max(model, closest_bdy_point);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);
        first_() = trace_ray_imaging<Rest>(model.geometry, origin, closest_bdy_point, ray_dir,
                       dshift_max, -1, Z, centre - 1, centre - 1)
                 + 1;
        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction

        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = image_shortchar_order_0<approx>(model, closest_bdy_point, f);
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = boundary_intensity(
                    model, closest_bdy_point, model.radiation.frequencies.nu(closest_bdy_point, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

///  Creates an image using the comoving solver (and the new image).
/// @param[in] model: the model to image
/// @param[in] ray_dir: the ray direction to image
/// @param[in] nxpix: number of pixels in the x-axis
/// @param[in] nypix: number of pixels in the y-axis
/// @param[in] image_freqs: list of frequencies to image
template <ApproximationType approx>
inline void Solver ::image_comoving_new_imager(Model& model, const Vector3D& ray_dir,
    const Size nxpix, const Size nypix, const Vector<Real>& image_freqs) {
    // note: we compute the comoving intensities using only a limited amount of frequencies and
    // interpolate the computed intensities afterwards

    Image image =
        Image(model.geometry, model.radiation.frequencies, Intensity, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);
    image.set_freqs(image_freqs); // when setting up the imager, it assumes by default that the
                                  // frequencies for the model are identical to the frequencies used
                                  // for the image. This is not the case for this comoving imager.

    setup_comoving_new_imager(model, image, ray_dir);

    std::multimap<Real, Size> template_multimap_image_freq_to_index;

    for (Size i; i < image_freqs.size(); i++) {
        template_multimap_image_freq_to_index.insert({image_freqs[i], i});
    }

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future;
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        // copy the image freq map, as the data will be consumed when interpolating the frequencies
        std::multimap<Real, Size> multimap_image_freq_to_index(
            template_multimap_image_freq_to_index);

        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);
        const Real dshift_max = get_dshift_max(model, closest_bdy_point);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);

        first_() = trace_ray_imaging<Rest>(model.geometry, origin, closest_bdy_point, ray_dir,
                       dshift_max, -1, Z, centre - 1, centre - 1)
                 + 1; // TODO: maybe instead compute the ray length, as we do not need to set
                      // anything here...

        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction
        n_tot_() = (last_() + 1) - first_();

        Vector<Size>& nr = nr_();

        Vector<double>& shift = shift_(); // contains the doppler shifts necessary for mapping
                                          // frequencies from observer to the comoving frame. As I
                                          // want to do the reverse, I need 2.0-shift instead.

        Size first_interesting_rayposidx = first_();
        Size last_interesting_rayposidx  = last_();

        // Now is the perfect time to setup the boundary conditions and data for the forward ray
        comoving_ray_bdy_setup_forward<approx>(model, last_interesting_rayposidx);

        Size rayposidx = first_() + 1; // ray position index -> point index through nr[rayposidx]
        // from first_+1 to the last interesting index, trace the ray
        while (rayposidx <= last_interesting_rayposidx) {
            const Real shift_next     = 2.0 - shift_()[rayposidx];
            const Real shift_curr     = 2.0 - shift_()[rayposidx - 1];
            const bool is_upward_disc = (shift_next >= shift_curr);

            solve_comoving_image_single_step(model, rayposidx, is_upward_disc, true);
            rayposidx++;
        }

        // We interpolate the computed intensities, starting from the last point
        rayposidx = last_();

        // Note that we can ignore the first point, as the intensities at this point are boundary
        // intensities
        while (rayposidx >= first_() + 1) {
            const Real shift_curr     = 2.0 - shift_()[rayposidx];
            const Real shift_prev     = 2.0 - shift_()[rayposidx - 1];
            const bool is_upward_disc = (shift_curr >= shift_prev);

            const Size currpoint = nr_()[rayposidx];

            // First, we figure out which line ranges exist for the current point, excluding the
            // boundary frequencies (interpolating them is not useful)
            get_line_ranges(model, currpoint, is_upward_disc, shift_curr);
            // And we check what boundary conditions need to be evaluated using the intensities at
            // curr_point
            interpolate_computed_comoving_intensities(model, image, pixidx, currpoint, rayposidx,
                shift_curr, multimap_image_freq_to_index);
            rayposidx--;
        }

        // All remaining frequencies have not been matched, so use boundary intensity for them
        for (std::pair<Real, Size> it : multimap_image_freq_to_index) {
            const Real image_freq = std::get<0>(it); // get image frequency
            const Size f          = std::get<1>(it); // get corresponding frequency index for image

            image.I(pixidx, f) = boundary_intensity(model, closest_bdy_point, image_freq);
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2_for_point(Model& model, const Size rr, const Size p) {
    // Redefine p to keep everything as similar to
    // image_feautrier_order_2 as possible
    const Size o = p;

    const Size ar = model.geometry.rays.get_antipod_index(rr);

    const Real dshift_max = get_dshift_max(model, o);

    nr_()[centre]    = o;
    shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);

    first_() =
        trace_ray<Rest, false>(model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1) + 1;
    last_() = trace_ray<Rest, false>(model.geometry, o, ar, dshift_max, +1, centre + 1, centre) - 1;
    n_tot_() = (last_() + 1) - first_();

    model.S_ray.resize(n_tot_(), model.parameters->nfreqs());
    model.dtau_ray.resize(n_tot_(), model.parameters->nfreqs());
    model.u_ray.resize(n_tot_(), model.parameters->nfreqs());

    if (n_tot_() > 1) {
        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            image_feautrier_order_2_for_point_loc<approx>(model, o, f);
        }
    }
}

template <ApproximationType approx>
inline void Solver ::image_optical_depth(Model& model, const Size rr) {
    Image image = Image(model.geometry, model.radiation.frequencies, Intensity, rr);

    accelerated_for(o, model.parameters->npoints(), {
        const Size ar = model.geometry.rays.get_antipod_index(rr);

        const Real dshift_max = get_dshift_max(model, o);

        nr_()[centre]    = o;
        shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);
        ;

        first_() =
            trace_ray<Rest, false>(model.geometry, o, rr, dshift_max, -1, centre - 1, centre - 1)
            + 1;
        last_() =
            trace_ray<Rest, false>(model.geometry, o, ar, dshift_max, +1, centre + 1, centre) - 1;
        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_optical_depth<approx>(model, o, f);

                image.I(o, f) = optical_depth_();
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(o, f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_optical_depth_new_imager(
    Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix) {
    Image image =
        Image(model.geometry, model.radiation.frequencies, OpticalDepth, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future; (but then while loop will be
    // required anyway)
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);
        const Real dshift_max = get_dshift_max(model, closest_bdy_point);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);
        first_() = trace_ray_imaging<Rest>(model.geometry, origin, closest_bdy_point, ray_dir,
                       dshift_max, -1, Z, centre - 1, centre - 1)
                 + 1;
        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction

        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_optical_depth<approx>(model, closest_bdy_point, f);

                image.I(pixidx, f) = optical_depth_();
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

// Because of the new method for computing the optical
// depth, adding extra frequency points for counteracting
// the large doppler shift is no longer necessary
template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::trace_ray(const Geometry& geometry, const Size o, const Size r,
    const double dshift_max, const int increment, Size id1, Size id2) {
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (geometry.valid_point(nxt)) {
        Size crt         = o;
        double shift_crt = geometry.get_shift<frame, use_adaptive_directions>(o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

        set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt)) {
            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            shift_nxt = geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

            set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

// Because of the new method for computing the optical
// depth, adding extra frequency points for counteracting
// the large doppler shift is no longer necessary
// Specialized raytracer for imaging the model. Note:
// assumes the outer boundary to be convex in order to
// correctly start tracing the ray
accel inline Size Solver ::trace_ray_imaging_get_start(const Geometry& geometry,
    const Vector3D& origin, const Size start_bdy, const Vector3D& raydir, Real& Z) {

    Size initial_point = start_bdy;
    // first figure out which boundary point lies closest to
    // the custom ray
    // TODO: is slightly inefficient implementation, can be
    // improved by only checking bdy point neighbors We use
    // here the assumption that the outer boundary is convex,
    // to obtain the closest point on the boundary
    while (true) {
        Size next_attempt =
            geometry.get_boundary_point_closer_to_custom_ray(origin, raydir, initial_point);
        if (next_attempt == initial_point) {
            break;
        }
        initial_point = next_attempt;
    }

    Z = geometry.get_distance_origin_to_boundary(origin, raydir, initial_point);

    return initial_point;
}

// Because of the new method for computing the optical
// depth, adding extra frequency points for counteracting
// the large doppler shift is no longer necessary
// Specialized raytracer for imaging the model. Note:
// assumes the outer boundary to be convex in order to
// correctly start tracing the ray
template <Frame frame>
accel inline Size Solver ::trace_ray_imaging(const Geometry& geometry, const Vector3D& origin,
    const Size start_bdy, const Vector3D& raydir, const double dshift_max, const int increment,
    Real& Z, // distance from origin can be non-zero to start,
             // as this measures the distance from the
             // projection plane
    Size id1, Size id2) {

    double dZ                      = 0.0; // last increment in Z
    const Vector3D origin_velocity = Vector3D(0.0, 0.0, 0.0);
    Size crt                       = start_bdy;

    Size nxt = geometry.get_next<Imagetracer>(origin, raydir, start_bdy, Z, dZ);

    if (geometry.valid_point(nxt)) {
        double shift_crt =
            geometry.get_shift<frame>(origin, origin_velocity, raydir, crt, Z, false);
        double shift_nxt =
            geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z, false);

        set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        // Due to boundary points begin slightly annoying to
        // start a ray from, we must make sure the ray only end
        // when we are sure that we stay on the boundary
        while (true) {
            if (!geometry.not_on_boundary(nxt)) {
                Size curr_cons_bdy = 1;
                Size temp_nxt      = nxt;
                double temp_Z      = Z;
                double temp_dZ     = dZ;

                while (true) {
                    temp_nxt =
                        geometry.get_next<Imagetracer>(origin, raydir, temp_nxt, temp_Z, temp_dZ);
                    if ((!geometry.valid_point(temp_nxt))
                        || (curr_cons_bdy == MAX_CONSECUTIVE_BDY)) {
                        return id1; // the ray ends if we cannot find
                                    // any points anymore, or we have
                                    // too many boundary points after
                                    // eachother
                    }
                    if (geometry.not_on_boundary(temp_nxt)) {
                        break; // the ray continues, as we find once
                               // again non-boundary points
                    }
                    curr_cons_bdy += 1;
                }
            }

            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = geometry.get_next<Imagetracer>(origin, raydir, nxt, Z, dZ);
            shift_nxt = geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z, false);

            set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

// Because of the new method for computing the optical depth, adding extra frequency points for
// counteracting the large doppler shift is no longer necessary

/// Traces the ray for the comoving solver, also returning the index of the outermost interesting
/// point
/// @param[in] geometry: model geometry to trace ray through
/// @param[in] o: point index of the starting point
/// @param[in] r: ray direction index for tracing the ray
/// @param[in] rr: ray direction index for checking whether the current ray lies closest to the
/// point in this raydirection; rr \in [0,geometry.parameters.hnrays()-1]
/// @param[in] rayidx: index of ray traced in this direction
/// @param[in] dshift_max: maximum allowed doppler shift (currently does not do anything)
/// @param[in] increment: increment for in which index to save the data
/// @param[in,out] id1: index at which to save traced point, and corresponding doppler shift
/// @param[in,out] id2: index at which to save distance increment (due to antipodal rays)
/// @param[out] outermost_interesting_point_rayidx: position index on ray of the last point for
/// which the traced ray lies closest to the point in the given ray direction
/// @return ray position index of last point put on ray + increment
/// TODO: deprecate rayidx usage
template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::trace_ray_comoving(const Geometry& geometry, const Size o, const Size r,
    const Size rr, const Size rayidx, const double dshift_max, const int increment, Size id1,
    Size id2, Size& outermost_interesting_point_rayidx) {
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (geometry.valid_point(nxt)) {
        Size crt         = o;
        double shift_crt = geometry.get_shift<frame, use_adaptive_directions>(o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

        // if (closest_ray(rr, nxt) == rayidx) {
        //     outermost_interesting_point_rayidx =
        //         id1; // as the data is set there, this should? be fine
        // }
        // Check whether we need to use this ray to compute some intensity data at point nxt
        if (intensity_origin[nxt].count(std::tuple(o, rr)))
        {
            outermost_interesting_point_rayidx = id1; // as the data is set there, this should? be fine
        }

        set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt)) {
            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            shift_nxt = geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

            // if (closest_ray(rr, nxt) == rayidx) {
            //     outermost_interesting_point_rayidx =
            //         id1; // as the data is set there, this should? be fine
            // }
            if (intensity_origin[nxt].count(std::tuple(o, rr)))
            {
                outermost_interesting_point_rayidx = id1; // as the data is set there, this should? be fine
            }

            set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

// Tracing the ray, ignoring all points for which the given CMF frequency lies too far from all line
// centers; TODO: either remove, or make compatible with imaging (different ray length for each
// frequency?)
template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::trace_ray_pruned(const Model& model, const Size o, const Size r,
    const double dshift_max, const int increment, Size id1, Size id2, const Real freq) {
    double Z     = 0.0; // distance from origin (o)
    double dZ    = 0.0; // last increment in Z
    double dZcrt = dZ;  // previous increment in dZ

    bool is_crt_set = true; // Whether the current point is already set

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt = o;
        // Size         prv = crt;
        double shift_crt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, crt, 0.0);
        // double shift_prv = shift_crt;
        double shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

        // As we must make sure that all lines are traced fine
        // check whether position increment has any close lines
        if (check_close_line(shift_crt * freq, shift_nxt * freq, crt, nxt, model))
        // if (check_close_line(shift_prv*freq, shift_crt*freq, shift_nxt*freq, prv, crt, nxt,
        // model))
        {
            set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
            is_crt_set = true;
        } else {
            is_crt_set = false;
        }

        while (model.geometry.not_on_boundary(nxt)) {
            //       prv =       crt;
            // shift_prv = shift_crt;
            crt       = nxt;
            shift_crt = shift_nxt;

            dZcrt = dZ; // current dZ required

            nxt       = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

            // As we must make sure that all lines are traced fine
            // check whether position increment has any close lines
            if (check_close_line(shift_crt * freq, shift_nxt * freq, crt, nxt, model))
            // if (check_close_line(shift_prv*freq, shift_crt*freq, shift_nxt*freq, prv, crt, nxt,
            // model))
            {
                if (!is_crt_set) {
                    // FIXME: FIX FUNCTION DEFINITION SET_DATA
                    set_data(
                        crt, crt, shift_crt, shift_crt, dZcrt, dshift_max, increment, id1, id2);
                }
                set_data(crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
                is_crt_set = true;
            } else {
                is_crt_set = false;
            }
        }

        // always set the boundary point (otherwise the boundary condition cannot be evaluated)
        if (!is_crt_set) {
            set_data(crt, crt, shift_crt, shift_crt, dZcrt, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

accel inline bool Solver ::check_close_line(const Real currfreq, const Real nextfreq,
    const Size currpoint, const Size nextpoint, const Model& model) {
    const Real left_freq  = std::min({currfreq, nextfreq});
    const Real right_freq = std::max({currfreq, nextfreq});

    // FIXME?: use wider bounds, as this is just for quantifying whether the line center lies close
    // enough. mhh, the farther we lie from the line center, the less influence the evaluated
    // frequency has on the mean line intensity... using maximum of bounds on the two points to get
    // an upper bound for the line width
    //  const Real prev_bound_line_width =
    //  model.thermodynamics.profile_width_upper_bound_with_linefreq(prevpoint, right_freq,
    //  model.lines.max_inverse_mass);
    const Real curr_bound_line_width = model.thermodynamics.profile_width_upper_bound_with_linefreq(
        currpoint, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.thermodynamics.profile_width_upper_bound_with_linefreq(
        nextpoint, right_freq, model.lines.max_inverse_mass);
    const Real upper_bound_line_width = model.parameters->max_distance_opacity_contribution
                                      * std::max({curr_bound_line_width, next_bound_line_width});

    const Real left_freq_bound  = left_freq - upper_bound_line_width;
    const Real right_freq_bound = right_freq + upper_bound_line_width;

    // apply default search algorithms on the bounds, obtaining iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);
    // yes, I am comparing pointers; I just want to know whether at least one line lies in the
    // interval // see compute_S_dtau_line_integrated <CloseLines>
    // TODO: refactor such that both trace_ray_pruned and compute_S_dtau_line_integrated
    // <CloseLines> rely on the same code snippet (instead of duplicating this code)
    return (left_line_bound != right_line_bound);
}

/// Computes the ray length when using the new imager
/// @param[in] geometry: geometry of the model to image
/// @param[in] origin: 3D position of the start of the ray
/// @param[in] start_bdy: point index of the first boundary point encountered by the ray
/// @param[in] raydir: ray direction vector of the ray
/// @return ray length
accel inline Size Solver ::get_ray_length_new_imager(const Geometry& geometry,
    const Vector3D& origin, const Size start_bdy, const Vector3D& raydir) {
    Size l             = 0;   // ray length, which we need to compute
    double Z           = 0.0; // total travelled distance
    double dZ          = 0.0; // last increment in Z
    Size initial_point = trace_ray_imaging_get_start(geometry, origin, start_bdy, raydir, Z);
    Size crt           = initial_point;
    Size nxt           = geometry.get_next<Imagetracer>(origin, raydir, crt, Z, dZ);

    if (geometry.valid_point(nxt)) {
        l += 1;

        // Due to boundary points begin slightly annoying to
        // start a ray from, we must make sure the ray only end
        // when we are sure that we stay on the boundary
        while (true) {
            if (!geometry.not_on_boundary(nxt)) {
                Size curr_cons_bdy = 1;
                Size temp_nxt      = nxt;
                double temp_Z      = Z;
                double temp_dZ     = dZ;
                while (true) {
                    temp_nxt =
                        geometry.get_next<Imagetracer>(origin, raydir, temp_nxt, temp_Z, temp_dZ);
                    if ((!geometry.valid_point(temp_nxt))
                        || (curr_cons_bdy == MAX_CONSECUTIVE_BDY)) {
                        return l; // the ray ends if we cannot find any
                                  // points anymore, or we have too many
                                  // boundary points after eachother
                    }
                    if (geometry.not_on_boundary(temp_nxt)) {
                        break; // the ray continues, as we find once
                               // again non-boundary points
                    }
                    curr_cons_bdy += 1;
                }
            }

            crt = nxt;
            nxt = geometry.get_next<Imagetracer>(origin, raydir, nxt, Z, dZ);
            l += 1;
        }
    }
    return l;
}

///  Sets the doppler shift, point index and position increment for a single position increment on a
///  ray. Also increments the index accordingly.
///   @param[in] crt: point index of the current point (currently no longer used)
///   @param[in] nxt: point index of the next point
///   @param[in] shift_crt: doppler shift of the current point (currently no longer used)
///   @param[in] shift_nxt: doppler shift of the next point
///   @param[in] dZ_loc: distance increment
///   @param[in] dshift_max: maximum doppler shift (currently no longer used)
///   @param[in] increment: index increment after setting the data
///   @param[in,out] id1: index at which to save traced point, and corresponding doppler shift
///   @param[in,out] id2: index at which to save distance increment (due to antipodal rays)
///   @note In the past, it was possible to split up a position increment into multiple interpolated
///   pieces in case of high velocity gradients. This was cumbersome to maintain and has been
///   removed as the sobolev-like treatment of optical depths can account for high velocity
///   gradients.
accel inline void Solver ::set_data(const Size crt, const Size nxt, const double shift_crt,
    const double shift_nxt, const double dZ_loc, const double dshift_max, const int increment,
    Size& id1, Size& id2) {
    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs(dshift);

    nr[id1]    = nxt;
    shift[id1] = shift_nxt;
    dZ[id2]    = dZ_loc;

    id1 += increment;
    id2 += increment;
}

///  Gaussian line profile function
///    @param[in] width : profile width
///    @param[in] diff  : frequency difference with line
///    centre
///    @return profile function evaluated with this
///    frequency difference
////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::gaussian(const Real inverse_width, const Real diff) const {
    const Real sqrt_exp = inverse_width * diff;

    return inverse_width * INVERSE_SQRT_PI * expf(-sqrt_exp * sqrt_exp);
}

///  Planck function
///    @param[in] temp : temperature of the corresponding
///    black body
///    @param[in] freq : frequency at which to evaluate the
///    function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::planck(const Real temp, const Real freq) const {
    return TWO_HH_OVER_CC_SQUARED * (freq * freq * freq) / expm1f(HH_OVER_KB * freq / temp);
}

///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] p      : point index of the boundary point
///    @param[in] freq   : frequency at which to evaluate
///    boundary condition
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::boundary_intensity(
    const Model& model, const Size p, const Real freq) const {
    const Size bdy_id = model.geometry.boundary.point2boundary[p];

    switch (model.geometry.boundary.boundary_condition[bdy_id]) {
    case Zero:
        return 0.0;
    case Thermal:
        return planck(model.geometry.boundary.boundary_temperature[bdy_id], freq);
    default:
        return planck(T_CMB, freq);
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<None>(const Model& model, const Size p,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters->nlines(); l++) {
        const Real diff = freq - model.lines.line[l];
        const Real prof = gaussian(model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity(p, l);
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///  function uses only the nearby lines to save computation
///  time
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<CloseLines>(const Model& model, const Size p,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    const Real upper_bound_line_width =
        model.parameters->max_distance_opacity_contribution
        * model.thermodynamics.profile_width_upper_bound_with_linefreq(
            p, freq, model.lines.max_inverse_mass);
    const Real left_freq_bound  = freq - upper_bound_line_width;
    const Real right_freq_bound = freq + upper_bound_line_width;

    // Just using default search algorithms, obtaining
    // iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++) {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // mapping sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];
        // const Real diff = freq - model.lines.line[l];
        const Real diff = freq - *freq_sort_l; // should be equal to the
                                               // previous line of code
        const Real inv_width = model.lines.inverse_width(p, l);
        const Real prof      = gaussian(model.lines.inverse_width(p, l), diff);
        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity(p, l);
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///  in the "one line" approixmation
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  l     : line index corresponding to the
///    frequency
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<OneLine>(
    const Model& model, const Size p, const Size l, const Real freq, Real& eta, Real& chi) const {
    const Real diff = freq - model.lines.line[l];
    const Real prof = gaussian(model.lines.inverse_width(p, l), diff);

    eta = prof * model.lines.emissivity(p, l);
    chi = prof * model.lines.opacity(p, l) + model.parameters->min_opacity;
}

///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap(const Real x_crt, const Real x_nxt, const double dZ) {
    return half * (x_crt + x_nxt) * dZ;
}

template <ApproximationType approx, bool use_adaptive_directions>
accel inline void Solver ::solve_shortchar_order_0(Model& model, const Size o, const Size r) {
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& source_c = source_c_();
    Vector<Real>& source_n = source_n_();

    Vector<Real>& tau = tau_();

    double Z  = 0.0; // distance along ray
    double dZ = 0.0; // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);
    Real term_c, term_n, dtau;
    bool compute_curr_opacity, prev_compute_curr_opacity;
    const Real ray_weight = model.geometry.rays.get_weight<use_adaptive_directions>(o, r);

    if (model.geometry.valid_point(nxt)) {
        double shift_c = 1.0;
        double shift_n = model.geometry.get_shift<CoMoving, use_adaptive_directions>(o, r, nxt, Z);

        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            const Real freq = model.radiation.frequencies.nu(o, f);
            const Size l    = model.radiation.frequencies.corresponding_line[f]; // line index

            compute_curr_opacity = true; // for the first point, we need to compute
                                         // both the curr and next opacity (and source)

            compute_source_dtau<approx>(model, crt, nxt, l, freq * shift_c, freq * shift_n, shift_c,
                shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f],
                source_n[f]);
            dtau = std::max(model.parameters->min_dtau, dtau);

            // proper implementation of 2nd order shortchar (not
            // yet times reducing factor of exp(-tau))
            //  model.radiation.I(r,o,f) = term_c *
            //  (expm1(-dtau)+dtau) / dtau
            //                           + term_n *
            //                           (-expm1(-dtau)-dtau*expf(-dtau))
            //                           /dtau;
            // Rewrite, trying to use less exponentials
            const Real factor = expm1f(-dtau) / dtau;

            model.radiation.I(r, o, f) =
                factor * (source_c[f] - source_n[f] * (1.0 + dtau)) + source_c[f] - source_n[f];
            tau[f] = dtau;

            // Compute local lambda operator
            const Size l_spec =
                model.radiation.frequencies.corresponding_l_for_spec[f]; // index of species
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[f]; // index of
                                                                                    // transition
            const Size z =
                model.radiation.frequencies.corresponding_z_for_line[f]; // index of
                                                                         // quadrature point
            const Real w_ang = ray_weight;

            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l_spec];

            const Real freq_line = lspec.linedata.frequency[k];
            const Real invr_mass = lspec.linedata.inverse_mass;
            const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;

            Real eta, chi; // eta is dummy var
            // chi is not necessarily computed, so compute it to
            // be sure
            get_eta_and_chi<approx>(model, o, k, freq_line, eta, chi);
            Real inverse_chi = 1.0 / chi;
            Real phi         = model.thermodynamics.profile(invr_mass, o, freq_line, freq);
            // const Real lambda_factor =
            // (dtau+expm1f(-dtau))/dtau;// If one wants to
            // compute lambda a bit more accurately in case of
            // dtau≃0. Real L   = constante * freq * phi *
            // lambda_factor * inverse_chi;
            Real L = constante * freq * phi * (factor + 1.0)
                   * inverse_chi; // using factor+1.0, the
                                  // computed lambda elements can
                                  // be negative if dtau very
                                  // small; but then the lambda
                                  // elements are also negligible
            lspec.lambda.add_element(o, k, o, L);

            // TODO: possible nonlocal lambda part // FIXME:
            // probably incorrect chi used
            //  L   = constante * freq * phi * (-factor *
            //  (1.0+dtau) - 1.0) * inverse_chi;
            //  lspec.lambda.add_element(o, k, nxt, L);
        }

        // For all frequencies, we need to use the same method
        // for computing the optical depth
        //  bool
        //  prev_compute_curr_opacity=compute_curr_opacity;//technically,
        //  we could also keep this bool individually for every
        //  frequency
        prev_compute_curr_opacity = compute_curr_opacity; // technically, we could also
                                                          // keep this bool individually
                                                          // for every frequency

        while (model.geometry.not_on_boundary(nxt)) {
            crt     = nxt;
            shift_c = shift_n;

            model.geometry.get_next<use_adaptive_directions>(o, r, crt, nxt, Z, dZ, shift_n);

            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                source_c[f]     = source_n[f];
                chi_c[f]        = chi_n[f];
                const Real freq = model.radiation.frequencies.nu(o, f);
                const Size l    = model.radiation.frequencies.corresponding_line[f];

                compute_curr_opacity = prev_compute_curr_opacity;

                compute_source_dtau<approx>(model, crt, nxt, l, freq * shift_c, freq * shift_n,
                    shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f],
                    source_c[f], source_n[f]);
                dtau = std::max(model.parameters->min_dtau, dtau);

                // proper implementation of 2nd order shortchar (not
                // yet times reducing factor of exp(-tau))
                //  model.radiation.I(r,o,f) += expf(-tau[f]) *
                //                           ( term_c *
                //                           (expm1(-dtau)+dtau) /
                //                           dtau
                //                           + term_n *
                //                           (-expm1(-dtau)-dtau*expf(-dtau))
                //                           /dtau);
                // Rewrite, trying to use less exponentials
                model.radiation.I(r, o, f) +=
                    expf(-tau[f])
                    * (expm1f(-dtau) / dtau * (source_c[f] - source_n[f] * (1.0 + dtau))
                        + source_c[f] - source_n[f]);
                // TODO: check order of addition, as we might be
                // starting with the largest contributions, before
                // adding the smaller ones...
                tau[f] += dtau;
            }

            // save setting for use for all frequencies for the
            // next interval
            prev_compute_curr_opacity = compute_curr_opacity;
        }

        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r, o, f) +=
                boundary_intensity(model, nxt, freq * shift_n) * expf(-tau[f]);
            model.radiation.J(o, f) += ray_weight * model.radiation.I(r, o, f);
        }
    }

    else {
        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r, o, f) = boundary_intensity(model, crt, freq);
            model.radiation.J(o, f) += ray_weight * model.radiation.I(r, o, f);
        }
    }
}

// Imager using long characteristics
// returns the computed specific intensity at the end of the ray
template <ApproximationType approx>
accel inline Real Solver ::image_shortchar_order_0(Model& model, const Size o, const Size f) {
    const Size l    = model.radiation.frequencies.corresponding_line[f];
    const Real freq = model.radiation.frequencies.nu(o, f);
    Real intensity  = 0.0;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();
    Size curr        = first;

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& source_c = source_c_();
    Vector<Real>& source_n = source_n_();

    Vector<Real>& tau = tau_();

    Size crt = nr[first];
    Size nxt = nr[first + 1];

    Real term_c, term_n, dtau;
    bool compute_curr_opacity, prev_compute_curr_opacity;

    if (last > curr) {
        double shift_c = shift[first];
        double shift_n = shift[first + 1];

        const Size l = model.radiation.frequencies.corresponding_line[f]; // line index

        compute_curr_opacity = true; // for the first point, we need to compute
                                     // both the curr and next opacity (and source)

        compute_source_dtau<approx>(model, crt, nxt, l, freq * shift_c, freq * shift_n, shift_c,
            shift_n, dZ[curr], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f],
            source_n[f]);
        dtau = std::max(model.parameters->min_dtau, dtau);

        // proper implementation of 2nd order shortchar (not
        // yet times reducing factor of exp(-tau))
        //  model.radiation.I(r,o,f) = term_c *
        //  (expm1(-dtau)+dtau) / dtau
        //                           + term_n *
        //                           (-expm1(-dtau)-dtau*expf(-dtau))
        //                           /dtau;
        // Rewrite, trying to use less exponentials
        const Real factor = expm1f(-dtau) / dtau;

        intensity = factor * (source_c[f] - source_n[f] * (1.0 + dtau)) + source_c[f] - source_n[f];
        tau[f]    = dtau;

        // For all frequencies, we need to use the same method
        // for computing the optical depth
        //  bool
        //  prev_compute_curr_opacity=compute_curr_opacity;//technically,
        //  we could also keep this bool individually for every
        //  frequency
        prev_compute_curr_opacity = compute_curr_opacity; // technically, we could also
                                                          // keep this bool individually
                                                          // for every frequency
        curr += 1;

        while (curr < last) {
            crt     = nr[curr];
            nxt     = nr[curr + 1];
            shift_c = shift[curr];
            shift_n = shift[curr + 1];

            // model.geometry.get_next(o, r, crt, nxt, Z, dZ, shift_n);

            // for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            source_c[f] = source_n[f];
            chi_c[f]    = chi_n[f];
            // const Real freq = model.radiation.frequencies.nu(o, f);
            // const Size l    = model.radiation.frequencies.corresponding_line[f];

            compute_curr_opacity = prev_compute_curr_opacity;

            compute_source_dtau<approx>(model, crt, nxt, l, freq * shift_c, freq * shift_n, shift_c,
                shift_n, dZ[curr], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f],
                source_n[f]);
            dtau = std::max(model.parameters->min_dtau, dtau);

            // proper implementation of 2nd order shortchar (not
            // yet times reducing factor of exp(-tau))
            //  model.radiation.I(r,o,f) += expf(-tau[f]) *
            //                           ( term_c *
            //                           (expm1(-dtau)+dtau) /
            //                           dtau
            //                           + term_n *
            //                           (-expm1(-dtau)-dtau*expf(-dtau))
            //                           /dtau);
            // Rewrite, trying to use less exponentials
            intensity += expf(-tau[f])
                       * (expm1f(-dtau) / dtau * (source_c[f] - source_n[f] * (1.0 + dtau))
                           + source_c[f] - source_n[f]);
            // TODO: check order of addition, as we might be
            // starting with the largest contributions, before
            // adding the smaller ones...
            tau[f] += dtau;
            // }

            // save setting for use for all frequencies for the
            // next interval
            prev_compute_curr_opacity = compute_curr_opacity;
            curr += 1;
        }

        intensity += boundary_intensity(model, nxt, freq * shift_n) * expf(-tau[f]);
        return intensity;
    }

    return boundary_intensity(model, crt, freq);
}

// accel inline void Solver ::
// solve_shortchar_order_0_ray_forward (
//           Model& model,
//           const Size   o,
//           const Size   r)
// {
//     Vector<Real>& eta_c = eta_c_();
//     Vector<Real>& eta_n = eta_n_();
//
//     Vector<Real>& chi_c = chi_c_();
//     Vector<Real>& chi_n = chi_n_();
//
//     Vector<Real>& source_c= source_c_();
//     Vector<Real>& source_n= source_n_();
//
//     Vector<Real>& tau = tau_();
//     Vector<double>& shift=shift_();
//     Vector<Size>& nr=nr_();
//     Vector<double>& dZ=dZ_();
//
//     Size crt, nxt;
//     // Size nxt = model.geometry.get_next (o, r, o, Z,
//     dZ); Real term_c, term_n, dtau; bool
//     compute_curr_opacity, prev_compute_curr_opacity;
//     prev_compute_curr_opacity=true;
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq =
//         model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model,
//         nr[first_()], freq*shift[first_()]);
//     }
//     double shift_c, shift_n;
//     //interate until we reach the middle point
//     for (Size idx=first_()+1; idx<=centre; idx++)
//     {
//
//         crt = nr[idx-1];
//         nxt = nr[idx];
//         shift_c = shift[idx-1];
//         shift_n = shift[idx];
//
//         for (Size f = 0; f < model.parameters->nfreqs();
//         f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq =
//             model.radiation.frequencies.nu(o, f); const
//             Size l    =
//             model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity =
//             prev_compute_curr_opacity; // for the first
//             point, we need to compute both the curr and
//             next opacity (and source)
//             // compute_curr_opacity = true; // for the
//             first point, we need to compute both the curr
//             and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l,
//             freq*shift_c, freq*shift_n, shift_c, shift_n,
//             dZ[idx-1], compute_curr_opacity, dtau,
//             chi_c[f], chi_n[f], source_c[f],
//             source_n[f]); dtau =
//             std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) =
//             exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c *
//             (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n *
//             (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) =
//             exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] *
//                                      (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] *
//                                        (expm1(-dtau)+dtau)
//                                        )/ dtau;
//         }
//         //save setting for use for all frequencies for
//         the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) +=
//         model.geometry.rays.weight[r] *
//         model.radiation.I(r,o,f);
//     }
// }
//
//
// accel inline void Solver ::
// solve_shortchar_order_0_ray_backward (
//           Model& model,
//           const Size   o,
//           const Size r)
// {
//     Vector<Real>& eta_c = eta_c_();
//     Vector<Real>& eta_n = eta_n_();
//
//     Vector<Real>& chi_c = chi_c_();
//     Vector<Real>& chi_n = chi_n_();
//
//     Vector<Real>& source_c= source_c_();
//     Vector<Real>& source_n= source_n_();
//
//     Vector<Real>& tau = tau_();
//     Vector<double>& shift=shift_();
//     Vector<Size>& nr=nr_();
//     Vector<double>& dZ=dZ_();
//
//     Size crt, nxt;
//     Real term_c, term_n, dtau;
//     bool compute_curr_opacity, prev_compute_curr_opacity;
//     prev_compute_curr_opacity=true; // for the first
//     point, we need to compute both the curr and next
//     opacity (and source)
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq =
//         model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model,
//         nr[last_()], freq*(shift[last_()]));
//     }
//     double shift_c, shift_n;
//     //interate until we reach the middle point
//     for (Size indexp1=last_(); indexp1>=centre+1;
//     indexp1--)
//     {
//         crt = nr[indexp1];
//         nxt = nr[indexp1-1];
//         shift_c = shift[indexp1];
//         shift_n = shift[indexp1-1];
//
//         for (Size f = 0; f < model.parameters->nfreqs();
//         f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq =
//             model.radiation.frequencies.nu(o, f); const
//             Size l    =
//             model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity =
//             prev_compute_curr_opacity; // for the first
//             point, we need to compute both the curr and
//             next opacity (and source)
//             // compute_curr_opacity = true; // for the
//             first point, we need to compute both the curr
//             and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l,
//             freq*shift_c, freq*shift_n, shift_c, shift_n,
//             dZ[indexp1-1], compute_curr_opacity, dtau,
//             chi_c[f], chi_n[f], source_c[f],
//             source_n[f]); dtau =
//             std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) =
//             exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c *
//             (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n *
//             (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) =
//             exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] *
//                                      (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] *
//                                        (expm1(-dtau)+dtau)
//                                        )/ dtau;
//         }
//
//         //save setting for use for all frequencies for
//         the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) +=
//         model.geometry.rays.weight[r] *
//         model.radiation.I(r,o,f);
//     }
// }

template <ApproximationType approx, bool use_adaptive_directions>
accel inline void Solver ::update_Lambda(Model& model, const Size rr, const Size f) {
    const Frequencies& freqs        = model.radiation.frequencies;
    const Thermodynamics& thermodyn = model.thermodynamics;

    if (freqs.appears_in_line_integral[f]) {
        const Size first = first_();
        const Size last  = last_();
        const Size n_tot = n_tot_();

        Vector<Size>& nr      = nr_();
        Vector<double>& shift = shift_();
        Vector<Real>& L_diag  = L_diag_();
        Matrix<Real>& L_upper = L_upper_();
        Matrix<Real>& L_lower = L_lower_();
        // Vector<Real  >& inverse_chi = inverse_chi_();

        const Real w_ang =
            two * model.geometry.rays.get_weight<use_adaptive_directions>(nr[centre], rr);

        const Size l = freqs.corresponding_l_for_spec[f]; // index of species
        const Size k = freqs.corresponding_k_for_tran[f]; // index of transition
        const Size z = freqs.corresponding_z_for_line[f]; // index of
                                                          // quadrature point

        LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

        const Real freq_line = lspec.linedata.frequency[k];
        const Real invr_mass = lspec.linedata.inverse_mass;
        const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;
        // TODO: approximation if we do not have overlapping
        // lines
        //  Real inverse_opacity = 1.0/model.lines.opacity
        //  (nr[centre], k);
        //  //Inverse line opacity; includes 1/HH_OVER_FOUR_PI
        //  Real L   = constante
        //  * L_diag[centre] * inverse_opacity;
        Real eta, chi; // eta is dummy var
        get_eta_and_chi<approx>(model, nr[centre], k, freq_line, eta, chi);
        Real inverse_chi = 1.0 / chi;

        Real frq = freqs.nu(nr[centre], f) * shift[centre];
        Real phi = thermodyn.profile(invr_mass, nr[centre], freq_line, frq);
        Real L   = constante * frq * phi * L_diag[centre] * inverse_chi;

        lspec.lambda.add_element(nr[centre], k, nr[centre], L);

        for (long m = 0; (m < n_off_diag) && (m + 1 < n_tot); m++) {
            if (centre >= first + m + 1) // centre-m-1 >= first
            {
                const long n = centre - m - 1;

                // TODO: approximation if we do not have overlapping
                // lines // also check which opacity we need
                // (centre?)
                //  inverse_opacity = 1.0/model.lines.opacity
                //  (nr[n], k); //Inverse line opacity; no longer
                //  1/HH_OVER_FOUR_PI L   = constante * L_lower(m,n)
                //  * inverse_opacity;
                get_eta_and_chi<approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi = 1.0 / chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile(invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_lower(m, n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }

            if (centre + m + 1 <= last) // centre+m+1 < last
            {
                const long n = centre + m + 1;

                // TODO: approximation if we do not have overlapping
                // lines
                //  inverse_opacity = 1.0/model.lines.opacity
                //  (nr[n], k); //Inverse line opacity; no includes
                //  1/HH_OVER_FOUR_PI L   = constante * L_upper(m,n)
                //  * inverse_opacity;
                get_eta_and_chi<approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi = 1.0 / chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile(invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_upper(m, n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }
        }
    }
}

///  Computer for the optical depth and source function when computing using the formal line
///  integration In case of low velocity differences, almost two times slower For high velocity
///  increments however, this does not need any extra interpolation points (-> way faster) (and some
///  extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its
///  limit value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<OneLine>(Model& model, Size currpoint,
    Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr,
    Real& Snext) {
    dtau  = compute_dtau_single_line(model, currpoint, nextpoint, lineidx, currfreq, nextfreq, dZ);
    Scurr = model.lines.emissivity(currpoint, lineidx)
          / model.lines.opacity(currpoint, lineidx); // current source
    Snext = model.lines.emissivity(nextpoint, lineidx)
          / model.lines.opacity(nextpoint, lineidx); // next source
    // note: due to interaction with dtau when computing all
    // sources individually, we do need to recompute Scurr and
    // Snext for all position increments
}

///  Computer for the optical depth and source function when computing using the formal line
///  integration In case of low velocity differences, almost two times slower For high velocity
///  increments however, this does not need any extra interpolation points (-> way faster) (and some
///  extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its
///  limit value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<None>(Model& model, Size currpoint,
    Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr,
    Real& Snext) {
    Real sum_dtau             = 0.0;
    Real sum_dtau_times_Scurr = 0.0;
    Real sum_dtau_times_Snext = 0.0;
    for (Size l = 0; l < model.parameters->nlines(); l++) {
        Real line_dtau =
            compute_dtau_single_line(model, currpoint, nextpoint, l, currfreq, nextfreq, dZ);
        Real line_Scurr = model.lines.emissivity(currpoint, l)
                        / model.lines.opacity(currpoint, l); // current source
        Real line_Snext =
            model.lines.emissivity(nextpoint, l) / model.lines.opacity(nextpoint, l); // next source
        sum_dtau += line_dtau;
        sum_dtau_times_Scurr += line_dtau * line_Scurr;
        sum_dtau_times_Snext += line_dtau * line_Snext;
    }
    dtau  = sum_dtau;
    Scurr = sum_dtau_times_Scurr / sum_dtau;
    Snext = sum_dtau_times_Snext / sum_dtau;
}

///  Computer for the optical depth and source function when
///  computing using the formal line integration In case of
///  low velocity differences, almost two times slower For
///  high velocity increments however, this does not need
///  any extra interpolation points (-> way faster) (and
///  some extra optimizations are made in the limit of
///  extremely large doppler shift (as erf goes to its limit
///  value) This function only takes into account the nearby
///  lines, saving some computation time
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<CloseLines>(Model& model, Size currpoint,
    Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr,
    Real& Snext) {
    Real sum_dtau = 0.0; // division by zero might occur otherwise
    // Real sum_dtau=model.parameters->min_dtau; //division by
    // zero might occur otherwise
    Real sum_dtau_times_Scurr = 0.0;
    Real sum_dtau_times_Snext = 0.0;

    Real left_freq;
    Real right_freq;

    // err, compiler will probably figure out that I just want
    // these two values ordered
    if (currfreq < nextfreq) {
        left_freq  = currfreq;
        right_freq = nextfreq;
    } else {
        right_freq = currfreq;
        left_freq  = nextfreq;
    }

    // using maximum of bounds on the two points to get an
    // upper bound for the line width
    const Real curr_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         currpoint, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         nextpoint, right_freq, model.lines.max_inverse_mass);
    const Real upper_bound_line_width = std::max(curr_bound_line_width, next_bound_line_width);

    const Real left_freq_bound  = left_freq - upper_bound_line_width;
    const Real right_freq_bound = right_freq + upper_bound_line_width;

    // apply default search algorithms on the bounds,
    // obtaining iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++) {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // Map sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];

        Real line_dtau =
            compute_dtau_single_line(model, currpoint, nextpoint, l, currfreq, nextfreq, dZ);
        Real line_Scurr = model.lines.emissivity(currpoint, l)
                        / model.lines.opacity(currpoint, l); // current source
        Real line_Snext =
            model.lines.emissivity(nextpoint, l) / model.lines.opacity(nextpoint, l); // next source
        sum_dtau += line_dtau;
        sum_dtau_times_Scurr += line_dtau * line_Scurr;
        sum_dtau_times_Snext += line_dtau * line_Snext;
    }
    dtau = sum_dtau;
    // needs extra bounding, as nothing may be added in the
    // first place (above for loop may have looped over 0
    // elements)
    const Real bound_min_dtau = model.parameters->min_opacity * dZ;
    // Correct way of bounding from below; should be able to
    // deal with very minor computation errors around 0.
    if (-bound_min_dtau < dtau) {
        dtau = std::max(bound_min_dtau, dtau);
    }
    // Note: 0 source functions can be returned if no lines
    // are nearby; but then the negligible lower bound gets
    // returned
    Scurr = sum_dtau_times_Scurr / dtau;
    Snext = sum_dtau_times_Snext / dtau;

    // note: due to interaction with dtau when computing all
    // sources individually, we do need to recompute Scurr and
    // Snext for all position increments
}

/// Computes the source function and optical depth in a
/// hybrid manner
///    @param[in/out] compute_curr_opacity: for deciding
///    whether we need to compute the current opacity when
///    using the trapezoidal rule
///    @param[in] currpoint : index of current point
///    @param[in] nextpoint : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] currshift : shift at curr point
///    @param[in] nextshift : shift at next point
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/// Warning: depending on how large the doppler shift is,
/// the opacity is NOT computed.
template <ApproximationType approx>
accel inline void Solver ::compute_source_dtau(Model& model, Size currpoint, Size nextpoint,
    Size line, Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ,
    bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr,
    Real& Snext) {
    // deciding which optical depth computation to use,
    // depending on the doppler shift
    const double dshift     = next_shift - curr_shift; // shift[first+1]-shift[first];TODO
    const double dshift_abs = fabs(dshift);
    const double dshift_max = std::min(model.dshift_max[currpoint], model.dshift_max[nextpoint]);
    const bool using_large_shift = (dshift_abs > dshift_max);

    // fancy computation for large doppler shifts
    if (using_large_shift) {
        compute_curr_opacity = true;
        compute_S_dtau_line_integrated<approx>(
            model, currpoint, nextpoint, line, curr_freq, next_freq, dZ, dtaunext, Scurr, Snext);
        // OPACITY IS NOT COMPUTED IN THIS BRANCH!
    } else {
        // default computation using trapezoidal rule
        if (compute_curr_opacity) // fancy computation does not
                                  // compute the current
                                  // opacity, so we might need
                                  // to recompute it here
        {
            compute_curr_opacity = false;
            Real eta_c           = 0.0; // current emissivity
            // also get previous opacity (emissivity does not
            // matter)
            get_eta_and_chi<approx>(model, currpoint, line, curr_freq, eta_c, chicurr);
            Scurr = eta_c / chicurr; // might as well compute the
                                     // source function too
        }

        Real eta_n = 0.0;
        // Get new radiative properties
        get_eta_and_chi<approx>(model, nextpoint, line, next_freq, eta_n, chinext);

        Snext = eta_n / chinext;

        dtaunext = half * (chicurr + chinext) * dZ;
    }
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::solve_feautrier_order_2(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], l, freq * shift[first],
        freq * shift[first + 1], shift[first], shift[first + 1], dZ[first], compute_curr_opacity,
        dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first]         = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = 1.0 / C[first]; // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], l, freq * shift[n],
            freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], compute_curr_opacity, dtau_n, chi_c,
            chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    if (n_off_diag == 0) {
        if (centre < last) {
            /// Write economically: G[last] = (B[last] - A[last])
            /// / A[last];
            GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
            GP[last] = GG[last] / (one + GG[last]);

            for (long n = last - 1; n > centre; n--) // use long in reverse loops!
            {
                Su[n] += Su[n + 1] * FI[n];

                GG[n] = (C[n] * GP[n + 1] + one) * inverse_A[n];
                GP[n] = GG[n] / (one + GG[n]);
            }

            Su[centre] += Su[centre + 1] * FI[centre];
            L_diag[centre] = inverse_C[centre] / (FF[centre] + GP[centre + 1]);
        } else {
            L_diag[centre] = (one + FF[centre - 1]) / (Bl_min_Al + Bl * FF[centre - 1]);
        }
    } else {
        /// Write economically: G[last] = (B[last] - A[last]) /
        /// A[last];
        GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
        GI[last] = one / (one + GG[last]);
        GP[last] = GG[last] * GI[last];

        L_diag[last] = (one + FF[last - 1]) / (Bl_min_Al + Bl * FF[last - 1]);

        for (long n = last - 1; n > first; n--) // use long in reverse loops!
        {
            Su[n] += Su[n + 1] * FI[n];

            GG[n] = (C[n] * GP[n + 1] + one) * inverse_A[n];
            GI[n] = one / (one + GG[n]);
            GP[n] = GG[n] * GI[n];

            L_diag[n] = inverse_C[n] / (FF[n] + GP[n + 1]);
        }

        Su[first] += Su[first + 1] * FI[first];
        L_diag[first] = (one + GG[first + 1]) / (Bf_min_Cf + Bf * GG[first + 1]);

        for (long n = last - 1; n >= first; n--) // use long in reverse loops!
        {
            L_upper(0, n + 1) = L_diag[n + 1] * FI[n];
            L_lower(0, n)     = L_diag[n] * GI[n + 1];
        }

        for (Size m = 1; (m < n_off_diag) && (m < n_tot - 1); m++) {
            for (long n = last - 1 - m; n >= first; n--) // use long in reverse loops!
            {
                L_upper(m, n + m + 1) = L_upper(m - 1, n + m + 1) * FI[n];
                L_lower(m, n)         = L_lower(m - 1, n) * GI[n + m + 1];
            }
        }
    }
}

///   Computes the optical depth assuming only a single line
///   exists.
///    @param[in] curridx: index of the current point
///    @param[in] nextidx: index of the next point
///    @param[in] lineidx: index of the line for which to
///    compute the optical depth
///    @param[in] curr_freq: current frequency (in comoving
///    frame)
///    @param[in] next_freq: next frequency (in comoving
///    frame)
///    @param[in] dz: distance increment
inline Real Solver ::compute_dtau_single_line(Model& model, Size curridx, Size nextidx,
    Size lineidx, Real curr_freq, Real next_freq, Real dz) {
    const Real linefreq = model.lines.line[lineidx];
    const Real average_inverse_line_width =
        (model.lines.inverse_width(curridx, lineidx) + model.lines.inverse_width(nextidx, lineidx))
        / 2.0;

    // opacity is stored divided by the linefreq, so multiply
    // by it
    const Real curr_line_opacity = model.lines.opacity(curridx, lineidx);
    const Real next_line_opacity = model.lines.opacity(nextidx, lineidx);

    // if frequencies are equal, division by zero (due to the
    // optical depth formula) happens if we were not to use
    // this branch
    if (curr_freq == next_freq) {
        // doing the default computation instead (no shifting)
        const Real diff = curr_freq - model.lines.line[lineidx]; // curr_freq==next_freq,
                                                                 // so choice is arbitrary
        const Real prof            = gaussian(average_inverse_line_width, diff);
        const Real average_opacity = (curr_line_opacity + next_line_opacity) / 2.0;

        return dz * (prof * average_opacity + model.parameters->min_opacity);
    }

    // We assume a linear interpolation of these dimensionless
    // frequency positions We will also assume the line width
    // to be somewhat constant, replacing the values with the
    // averages
    const Real next_pos = (linefreq - next_freq) * average_inverse_line_width;
    const Real curr_pos = (linefreq - curr_freq) * average_inverse_line_width;

    // In this way, the diff_pos can be computed quite simple,
    // and we do not have a discrepancy between the
    // interpolation and the bounds
    const Real diff_pos = next_pos - curr_pos;

    /// the more correct approach, taking into account also
    /// the line opacity change; however, it does not make too
    /// much of a difference in the actual result and is quite
    /// a bit slower
    // const Real
    // delta_opacity=(next_line_opacity-curr_line_opacity);
    // const Real deltanu=-next_freq+curr_freq;//differences
    // in curr freqs; +-1 due to shift being defined in the
    // other direction
    //
    // //note: opacity can also be extrapolated; however the
    // correction term (expterm) accounts for that const Real
    // interp_opacity=curr_line_opacity+delta_opacity*(curr_freq-linefreq)/deltanu;
    //
    // //This term is a constant term, giving the usual ... as
    // if the opacity were static const Real
    // erfterm=interp_opacity/diff_pos/2.0*(std::erf(next_pos)-std::erf(curr_pos));
    // //This term corrects for the fact that the opacity
    // between points changes const Real
    // expterm=delta_opacity/2.0*INVERSE_SQRT_PI/diff_pos/diff_pos*(std::exp(-curr_pos*curr_pos)-std::exp(-next_pos*next_pos));
    // return
    // dz*std::max(average_inverse_line_width*(erfterm+expterm),
    // model.parameters->min_opacity);

    // If we instead use an average opacity, the computation
    // is quite a bit faster
    const Real average_opacity = (next_line_opacity + curr_line_opacity) / 2.0;
    const Real erfterm =
        average_opacity / diff_pos / 2.0 * (std::erff(next_pos) - std::erff(curr_pos));
    // correcting to bound opacity from below to the minimum
    // opacity (assumes positive opacities occuring in the
    // model)
    return dz * (average_inverse_line_width * erfterm + model.parameters->min_opacity);
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_feautrier_order_2(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], l, freq * shift[first],
        freq * shift[first + 1], shift[first], shift[first + 1], dZ[first], compute_curr_opacity,
        dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], l, freq * shift[n],
            freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], compute_curr_opacity, dtau_n, chi_c,
            chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    // for (long n = last-1; n > first; n--) // use long in
    // reverse loops!
    // {
    //     Su[n] += Su[n+1] * FI[n];
    // }

    // Su[first] += Su[first+1] * FI[first];
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_feautrier_order_2_for_point_loc(
    Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], l, freq * shift[first],
        freq * shift[first + 1], shift[first], shift[first + 1], dZ[first], compute_curr_opacity,
        dtau_n, chi_c, chi_n, term_c, term_n);

    // err, source function might be slightly different when
    // looking at it from curr and next point
    //  this is due to the weighting by the line optical
    //  depths; this might not be saved TODO think whether
    //  this is correct
    model.S_ray(0, f) = term_c;
    model.S_ray(1, f) = term_n;

    model.dtau_ray(0, f) = 0.0;
    model.dtau_ray(1, f) = dtau_n;

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], l, freq * shift[n],
            freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], compute_curr_opacity, dtau_n, chi_c,
            chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        model.S_ray(n + 1 - first, f)    = term_n;
        model.dtau_ray(n + 1 - first, f) = dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    model.u_ray(last - first, f) = Su[last];

    for (long n = last - 1; n > first; n--) // use long in reverse loops!
    {
        Su[n] += Su[n + 1] * FI[n];

        model.u_ray(n - first, f) = Su[n];
    }

    Su[first] += Su[first + 1] * FI[first];

    model.u_ray(0, f) = Su[first];
}

/// Image the optical depth
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_optical_depth(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    Real tau = 0.0;

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], l, freq * shift[first],
        freq * shift[first + 1], shift[first], shift[first + 1], dZ[first], compute_curr_opacity,
        dtau_n, chi_c, chi_n, term_c, term_n);
    tau += dtau_n;

    /// iterate over all other points on the ray
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], l, freq * shift[n],
            freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], compute_curr_opacity, dtau_n, chi_c,
            chi_n, term_c, term_n);
        tau += dtau_n;
    }

    optical_depth_() = tau;
}

///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///  Computes both the mean intensity u and the flux v
//////////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::solve_feautrier_order_2_uv(Model& model, const Size o, const Size f) {
    // Note: in comments, we have another option for computing v (Sv at every location), by
    // mirroring the solution steps for u with slightly different data. This might be slower, and a
    // bit more work to make second-order accurate.
    // Currently, only sv at the center position gets computed using du/dtau = -v
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_();
    Vector<Size>& nr      = nr_();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();

    bool compute_curr_opacity =
        true; // for the first point, we need to compute both the curr and next opacity(and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], l, freq * shift[first],
        freq * shift[first + 1], shift[first], shift[first + 1], dZ[first], compute_curr_opacity,
        dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first]         = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = one / C[first]; // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);
    // TODO: if re-implementing other option, make this dSdtau second order accurate
    // Real dSdtau          = (term_n - term_c) / dtau_n;
    // only first order accurate dSdtau, as otherwise the solver needs to be rewritten
    // Current algorithm only allows us to access the data at current and next point

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    // Sv[first] = -dSdtau + two * inverse_dtau_f * (I_bdy_f - term_c);

    Su[first] /= Bf;
    // Sv[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], l, freq * shift[n],
            freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], compute_curr_opacity, dtau_n, chi_c,
            chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        // TODO: if re-implementing other option, make this dSdtau second order accurate
        // dSdtau              = (term_n - term_c) / dtau_n;
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
        // Sv[n] = (A[n] * Sv[n - 1] - dSdtau) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    // Sv[last] = -dSdtau - two * inverse_dtau_l * (I_bdy_l - term_n);
    // Different sign for Sv last boundary condition extra term (2/dtau(I-S))! (should be
    // assymetric)

    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;
    // Sv[last] = (A[last] * Sv[last - 1] + Sv[last]) * (one + FF[last - 1]) * denominator;

    if (centre < last) {
        for (long n = last - 1; n >= centre; n--) // use long in reverse loops !
        {
            Su[n] += Su[n + 1] * FI[n];
            // Sv[n] += Sv[n + 1] * FI[n];
        }
    } else {
        // Compute v using boundary at the end
        Sv[last] = Su[last] - I_bdy_l;
        return;
    }
    // Do one extra step for computing the derivative of the mean intensity
    if (centre > first) {
        Su[centre - 1] += Su[centre] * FI[centre - 1];
        // Recompute the optical depth increments, as we forgot to save them; FIXME: add
        // threadprivate dtau_() to solver.hpp
        compute_curr_opacity = true;
        compute_source_dtau<approx>(model, nr[centre - 1], nr[centre], l, freq * shift[centre - 1],
            freq * shift[centre], shift[centre - 1], shift[centre], dZ[centre - 1],
            compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
        const Real dtaumin = dtau_n;
        compute_source_dtau<approx>(model, nr[centre], nr[centre + 1], l, freq * shift[centre],
            freq * shift[centre + 1], shift[centre], shift[centre + 1], dZ[centre],
            compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
        const Real dtauplus = dtau_n;
        // TODO: optimize this calculation
        const Real coeffmin  = -dtauplus / (dtaumin * dtaumin + dtaumin * dtauplus);
        const Real coeffplus = dtaumin / (dtaumin * dtauplus + dtauplus * dtauplus);
        const Real coeffzero = -coeffmin - coeffplus;

        Sv[centre] =
            -(coeffmin * Su[centre - 1] + coeffplus * Su[centre + 1] + coeffzero * Su[centre]);
    } else {
        // Compute v using boundary at the start
        Sv[first] = I_bdy_f - Su[first];
        return;
    }
}

accel inline void Solver ::set_eta_and_chi(Model& model, const Size rr) const {
    model.eta.resize(model.parameters->npoints(), model.parameters->nfreqs());
    model.chi.resize(model.parameters->npoints(), model.parameters->nfreqs());

    for (Size p = 0; p < model.parameters->npoints(); p++) {
        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            // Extract the Doppler shift
            const double shift = model.geometry.get_shift<Rest, false>(0, rr, p, 0.0);
            const Real freq    = model.radiation.frequencies.nu(0, f);
            const Size l       = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi<None>(model, p, l, freq * shift, model.eta(p, f), model.chi(p, f));
        }
    }
}

accel inline void Solver ::set_boundary_condition(Model& model) const {
    model.boundary_condition.resize(model.parameters->nboundary(), model.parameters->nfreqs());

    for (Size b = 0; b < model.parameters->nboundary(); b++) {
        const Size p = model.geometry.boundary.boundary2point[b];

        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            const Real freq = model.radiation.frequencies.nu(0, f);

            model.boundary_condition(b, f) = boundary_intensity(model, p, freq);
        }
    }
}

template <bool use_adaptive_directions> inline void Solver ::set_column(Model& model) const {
    model.column.resize(model.parameters->nrays(), model.parameters->npoints());

    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar       = model.geometry.rays.get_antipod_index(rr);
            model.column(rr, o) = get_column<use_adaptive_directions>(model, o, rr);
            model.column(ar, o) = get_column<use_adaptive_directions>(model, o, ar);
        })
    }
}

template <bool use_adaptive_directions>
accel inline Real Solver ::get_column(const Model& model, const Size o, const Size r) const {
    Real column = 0.0;

    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt = o;

        column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;

        while (model.geometry.not_on_boundary(nxt)) {
            crt = nxt;
            nxt = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);

            column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;
        }
    }

    return column;
}
