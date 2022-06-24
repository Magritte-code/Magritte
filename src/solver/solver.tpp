template <Frame frame>
inline void Solver :: setup (Model& model)
{
    const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size  width = model.parameters->nfreqs();
    const Size  n_o_d = model.parameters->n_off_diag;

    setup (length, width, n_o_d);
}


inline void Solver :: setup (const Size l, const Size w, const Size n_o_d)
{
    length     = l;
    centre     = l/2;
    width      = w;
    n_off_diag = n_o_d;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        eta_c_       (i).resize (width);
        eta_n_       (i).resize (width);

        chi_c_       (i).resize (width);
        chi_n_       (i).resize (width);

        inverse_chi_ (i).resize (length);

        tau_         (i).resize (width);

        Su_          (i).resize (length);
        Sv_          (i).resize (length);

        A_           (i).resize (length);
        C_           (i).resize (length);
        inverse_A_   (i).resize (length);
        inverse_C_   (i).resize (length);

        FF_          (i).resize (length);
        FI_          (i).resize (length);
        GG_          (i).resize (length);
        GI_          (i).resize (length);
        GP_          (i).resize (length);

        L_diag_      (i).resize (length);

        L_upper_     (i).resize (n_off_diag, length);
        L_lower_     (i).resize (n_off_diag, length);
    }
}


// /  Getter for the maximum allowed shift value determined by the smallest line
// /    @param[in] o : number of point under consideration
// /    @retrun maximum allowed shift value determined by the smallest line
// /////////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: get_dshift_max (
    const Model& model,
    const Size   o     )
{
    Real dshift_max = std::numeric_limits<Real>::max();

    for (const LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        const Real inverse_mass   = lspec.linedata.inverse_mass;
        const Real new_dshift_max = model.parameters->max_width_fraction
                                    * model.thermodynamics.profile_width (inverse_mass, o);

        if (dshift_max > new_dshift_max)
        {
            dshift_max = new_dshift_max;
        }
    }

    return dshift_max;
}


template <Frame frame>
inline void Solver :: get_ray_lengths (Model& model)
{
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            model.geometry.lengths(rr,o) =
                model.geometry.get_ray_length <frame> (o, rr, dshift_max)
              + model.geometry.get_ray_length <frame> (o, ar, dshift_max);
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec();
}


template <Frame frame>
inline Size Solver :: get_ray_lengths_max (Model& model)
{
    get_ray_lengths <frame> (model);

    Geometry& geo = model.geometry;

    geo.lengths_max = *std::max_element(geo.lengths.vec.begin(),
                                        geo.lengths.vec.end()   );

    return geo.lengths_max;
}

//traces all points on a given ray (should be called in both directions)
//but as tracing antipodal rays should be symmetric, we only need to keep track of one direction
accel inline void Solver :: trace_ray_points (
    const Geometry& geometry,
    const Size      o,
    const Size      rdir,
    const Size      rsav,
    const Size      rayidx)
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, rdir, o, Z, dZ);

    // // For generality, I assign a ray index to each ray of a given direction rsav
    // // TODO: maybe change this if we were to allow arbitrary rays to be traced
    // // Also the ray direction is identified with the forward dir
    // closest_ray(rsav,o)=rayidx;
    // min_ray_distsqr(rsav,o)=0.0;

    if (geometry.valid_point(nxt))
    {
        n_rays_through_point(rsav,nxt)++;

        // get distance and check if closest ray
        // TODO: maybe write optimize function which also returns the distance, as we are currently doing some minor work twice
        Real dist2 = geometry.get_dist2_ray_point(o, nxt, r)

        if (dist2<min_ray_distsqr(rsav,nxt))
        {
            min_ray_distsqr(rsav,nxt)=dist2;
            closest_ray(rsav,nxt)=rayidx;
        }

        Size         crt = o;

        while (geometry.not_on_boundary(nxt))
        {
            crt =       nxt;
            nxt = geometry.get_next(o, rdir, nxt, Z, dZ);

            n_rays_through_point(rsav,nxt)++;

            // get distance and check if closest ray
            Real dist2 = geometry.get_dist2_ray_point(o, nxt, r)

            if (dist2<min_ray_distsqr(rsav,nxt))
            {
                min_ray_distsqr(rsav,nxt)=dist2;
                closest_ray(rsav,nxt)=rayidx;
            }
        }
    }
}

inline void Solver :: get_static_rays_to_trace (Model& model)
{
    accelerated_for (rr, model.parameters.hnrays(),
    {
        Size n_rays_to_trace=0;
        const Size ar=model.geometry.rays.antipod[rr];
        for (Size o=0; o<model.parameters.npoints(); o++)
        {
            // std::cout<<"n_rays_through_point: "<<n_rays_through_point(rr,o)<<std::endl;
            // TODO?: maybe replace with some boolean-like thing, as we do not actually care (aside from debugging)
            // how much rays are traced through a point
            if (n_rays_through_point(rr,o)>0)
            {continue;}

            //seemingly no ray has been traced through this point, so we must trace a ray through it
            points_to_trace_ray_through[rr][n_rays_to_trace]=o;
            n_rays_to_trace++;
            //tracing rays is symmetric, so only keep for the first half of the ray directions

            //trace ray through point
            n_rays_through_point(rr,o)++;
            //antipod has exactly same number of rays through the point, so do not save

            // For generality, I assign a ray index to each ray of a given direction rsav
            // Also the ray direction is identified with the forward dir
            closest_ray(rsav,o)=rayidx;
            min_ray_distsqr(rsav,o)=0.0;

            //now trace rest of rays
            trace_ray_points(model.geometry, o, rr, rr);
            trace_ray_points(model.geometry, o, ar, rr);
        }
        // n_points_to_trace_ray_through[rr]=n_rays_to_trace;
        points_to_trace_ray_through[rr].resize(n_rays_to_trace);//and now the correct size, instead of parameters->npoints()

    })

    //debug print stuff
    for (Size rr=0; rr<model.parameters.hnrays(); rr++)
    {
        std::cout<<"rr: "<<rr<<" size points_to_trace_ray_through: "<<points_to_trace_ray_through[rr].size()<<std::endl;
        for (Size idx=0; idx<points_to_trace_ray_through[rr].size(); idx++)
        {
            // std::cout<<"point: "<<points_to_trace_ray_through[rr][idx]<<std::endl;
        }
        // std::cout<<"number of rays per point"<<std::endl;
        for (Size p=0; p<model.parameters.npoints(); p++)
        {
            // std::cout<<"point: "<<p<<"#: "<<n_rays_through_point(rr, p)<<std::endl;
        }

    }

}


// Splits the frequencies into different parts, such that the different parts are sufficiently seperated
// TODO: In single line approx, this step should be ignored, as all lines are far enough from eachother
//ASSUMES FREQS SORTED
// p should be next point on the ray
//NOTE TO SELF: I ASSUMED FOR SIMPLICITY THAT I SPLIT IT A SINGLE RAY AT A TIME; NOT DIFFERENTLY FOR EACH POINT ON THE RAY
//TODO FIXME: USE THE RAY TO SPLIT THE FREQUENCIES IN A GLOBAL BOUND MANNER
// as a consequence, handling frequency quadratures very close to eachother should be automatically included
//Note: this can mean however that the frequency derivative terms might not be computed that accurately if the lines separate
inline void Solver :: split_frequencies (Model& model, Size p)
{
    freqsplits_()[0]=0;
    n_freqsplits_()=1;
    Real currfreq=model.radiation.frequencies.nu(p, 0);
    Real deltafreq=0.0;
    //TODO GET MIN LINE WIDTH
    for (freqid=1; freqid<model.parameters->nfreqs(), freqid++)
    {
        deltafreq=model.radiation.frequencies.nu(p, freqid)-currfreq;
        if (deltafreq>FACTORTODO*LINEWIDTHTODO)
        {
            freqsplits_()[n_freqsplits_()]=freqid;
            n_freqsplits_()++;
        }
        currfreq=model.radiation.frequencies.nu(p, freqid);
    }
    freqsplits_()[n_freqsplits_()]=model.parameters->nfreqs();//for simplicity, also add the last split.. makes it easier to get rightmost value in a split
}


inline void Solver :: solve_shortchar_order_0 (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.I.resize (model.parameters->nrays(),  model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize (                           model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Initialise mean intensity
    model.radiation.initialize_J();


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            // const Real dshift_max = get_dshift_max (o);
            const Real dshift_max = 1.0e+99;

            solve_shortchar_order_0 (model, o, rr, dshift_max);
            solve_shortchar_order_0 (model, o, ar, dshift_max);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                // model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_uv (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    solve_feautrier_order_2_uv <approx> (model, o, f);

                    model.radiation.u(rr,o,f)  = Su_()[centre];
                    model.radiation.v(rr,o,f)  = Sv_()[centre];
                }
            }
            else
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.v(rr,o,f)  = 0.0;
                }
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.v.copy_ptr_to_vec();
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_sparse (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr] * two;
        const Vector3D nn = model.geometry.rays.direction[rr];

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);

                nr_   ()[centre] = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                n_tot_() = (last_()+1) - first_();

                if (n_tot_() > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o,  lspec.nr_line[o][k][z]);

                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda           (model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}

// Matches the frequency indices such that the frequency difference is minimal
// ALSO relies on the current frequency split
// ERR, should just iterate over all frequencies, ignoring those different splits. Boundary conditions
// is_upward_disc denotes whether we want an upward discretization
// nextpointonrayindex should denote how far the next point lies on the ray
// currpointonrayindex should do the same, but for the previous point
inline void Solver :: match_frequency_indices(Model& model, const Size nextpoint, const Size currpoint, Size nextpointonrayindex, Size currpointonrayindex, bool is_upward_disc)
{//, const Size splitindex
    //if upward discretization, we start from the uppermost part
    if (is_upward_disc)
    {
        //Starting from the highest frequency
        //+1 to all indices due to using unsigned ints when looping down (overflow otherwise)
        Size curr_freq_idx=model.parameters->nfreqs()-1+1;
        Size next_freq_idx=model.parameters->nfreqs()-1+1;
        const Size min_freq_idx=0+1;
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        //NOTE: oob freqs should be dealt with later on
        while (curr_freq_idx>=min_freq_idx && next_freq_idx>=min_freq_idx)
        {
            //check if a the frequency at the previous point index is higher than this frequency
            if (model.radiation.frequencies.nu(currpoint, curr_freq_idx-1)>
                model.radiation.frequencies.nu(nextpoint, next_freq_idx-1))
            {
                curr_freq_idx--;
                // freq_matching_()[next_freq_idx]=curr_freq_idx;
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                // freq_matching_()[next_freq_idx-1]=curr_freq_idx-1;
                start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point on ray index
                start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=curr_freq_idx-1;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx--;
            }
        }
        //match remaining points as best as possible (with the outermost point)
        while (next_freq_idx>=min_freq_idx)
        {
            start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point
            start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=min_freq_idx-1;//stores freq
            next_freq_idx--;
        }


    }
    else
    {
        //Starting from the lowests frequency
        Size curr_freq_idx = 0;
        Size next_freq_idx = 0;
        const Size max_freq_idx = model.parameters->nfreqs()-1;
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        while (curr_freq_idx<=max_freq_idx && temp_freq_idx<=max_freq_idx)
        {
            //check if a previous index exists which is higher than this
            if (model.radiation.frequencies.nu(currpoint, curr_freq_idx)<
                model.radiation.frequencies.nu(nextpoint, next_freq_idx))
            {
                curr_freq_idx++;
                // freq_matching_()[next_freq_idx]=curr_freq_idx;
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                // freq_matching_()[next_freq_idx]=curr_freq_idx;
                start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
                start_indices_()(nextpointonrayindex,next_freq_idx)[1]=curr_freq_idx;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx++;
            }
        }
        //match remaining points as best as possible (with the outermost point)
        while (next_freq_idx<=max_freq_idx)
        {
            start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
            start_indices_()(nextpointonrayindex,next_freq_idx)[1]=max_freq_idx;//stores freq
            next_freq_idx++;
        }
    }
//? should also return the number/indices of boundary frequencies
}

// //TODO: not really needed anymore, as the set_[left/right]_implicit_boundary_frequencies does the same job, but better suited
// //TODO? seperate the curr and next boundary conditions? No, just do it in between the explicit and implicit part
// inline void Solver :: set_boundary_frequencies(const Size nextpoint, const Size currpoint, const Size splitindex, bool is_upward_disc)
// {
//     //TODO FOR SPLITS
//     //somehow intersect the frequency ranges of the next point and previous point, given a given split.
//     //split min, max indices are given by
//     const Size min_freq_idx=freqsplits_()[splitindex];//freq index for the current point
//     const Size max_freq_idx = freqsplits_()[splitindex+1]-1//minimal freq index for this split
//
//     const Real curr_min_bound = model.radiation.frequencies.nu(currpoint, min_freq_idx);
//     const Real curr_max_bound = model.radiation.frequencies.nu(currpoint, max_freq_idx);
//     const Real next_min_bound = model.radiation.frequencies.nu(nextpoint, min_freq_idx);
//     const Real next_max_bound = model.radiation.frequencies.nu(nextpoint, max_freq_idx);
//
//     //lazy checking whether ourside the bounds; TODO improve, as this currently does not use the fact that freqs are sorted
//     for (Size curr_freq_idx=min_freq_idx; curr_freq_idx<=max_freq_idx; curr_freq_idx++)
//     {
//         const Real temp_curr_freq = model.radiation.frequencies.nu(currpoint, curr_freq_idx);
//         if (temp_curr_freq<next_min_bound || temp_curr_freq>next_max_bound)
//         {
//             //TODO DATASTRUCT OR JUST DIRECTLY DO THE MOVE TO BDY INTENSITIES
//
//         }
//     }
//
//     //lazy checking whether ourside the bounds; TODO improve, as this currently does not use the fact that freqs are sorted
//     for (Size next_freq_idx=min_freq_idx; next_freq_idx<=max_freq_idx; next_freq_idx++)
//     {
//         const Real temp_next_freq = model.radiation.frequencies.nu(nextpoint, next_freq_idx);
//         if (temp_next_freq<curr_min_bound || temp_next_freq>curr_max_bound)
//         {
//             //TODO DATASTRUCT OR JUST DIRECTLY GET FROM BDY INTENSITIES
//             is_bdy_freq_()[temp_next_freq]=true;
//             //Also set boundary frequency
//             //using curr_min_bound or curr_max_bound as values for interpolating
//             //remove from list when encountering smaller/larger / ALSO EQUAL values as stored
//         }
//     }
//
//     //TODO
//
//     //for next point, set boundary points to background/prev intensity
//     //for curr point, store boundary points
// }

//Computes which lines are overlapping at the given point (next point)
//might not be useful, better to determine covering; although we might need it anyway for determining the freq der boundary conditions
//TODO: determine both (and define data struct for the overlap stuff) technically a simple Vector should work fine (even/odd indices respectively denotes start and end of intervals)
//err, we do need them at different points (overlap at curr point, freq der at next point hmm, also at curr?)
//USING THE SINGLE LINE APPROX, THIS CAN BE SIMPLIFIED
inline void Solver :: get_overlapping_lines(Model& model, const Size nextpoint, bool is_upward_disc)
{
    //Copied logic from get_line_ranges, as this should do similar things
    //Should be conceptually more simple than using the line centers and widths to figure out whether they overlap
    Vector<unsigned char>& line_quad_discdir_overlap=line_quad_left_overlap_();
    for (Size lineidx=0; lineidx<model.parameters->nlines(); lineidx++)
    {
        line_quad_discdir_overlap[lineidx]=false;//initialize to no overlap
    }

    Vector<Size>& line_count=line_count_();//counts the number of quadratures encountered for each line
    //This can get complicated due to having no requirement to use the same quadrature for all lines...
    //Due to having no requirement that all lines should use the same quadrature, we will just try to compute continguous ranges by iterating over the frequencies itself.
    //Thus we iterate over all freqs, and count the occurences of the lines; if we encounter a line the maximal number of times, set line count to zero and check if continguous region has ended
    // if not, then we evidently have an overlap
    //This only works because the frequencies themselves are sorted
    if (is_upward_disc)
    {
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            //get corresponding line of the freq
            const Size lineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, freqidx);
            //increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            //if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx]==model.parameters->nfreqs())
            {
                line_count[lineidx]=0;//reset the index (as counting this line is no longer needed)
                //Then if no other lines are currently being counted (i.e. entire count vectors is zero),
                // we conclude that the current continguous range has ended
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
                {
                    //just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx]=true;
                }
            }
        }
    }
    else
    {   //long for more easily handling reverse loops
        for (long freqidx=model.parameters->nfreqs()-1; freqidx>=0; freqidx--)
        {
            //get corresponding line of the freq
            const Size lineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, freqidx);
            //increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            //if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx]==model.parameters->nfreqs())
            {
                line_count[lineidx]=0;//reset the index (as counting this line is no longer needed)
                //Then if no other lines are currently being counted (i.e. entire count vectors is zero),
                // we conclude that the current continguous range has ended
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
                {
                    //just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx]=true;
                }
            }
        }
    }
    //TODO: check line_count to be uniformly 0!!
}

//Computes the continguous ranges spanned by the line quadratures
//excludes some part of the frequency quadrature (2 farthest freq quads) due to imposed boundary conditions
//USING THE SINGLE LINE APPROX, THIS CAN BE SIMPLIFIED ENORMOUSLY (just use all individual line quadratures)
inline void Solver :: get_line_ranges(Model& model, const Size curr_point, bool is_upward_disc, Real curr_shift)
{
    Vector<Real>& left_bound=left_bound_();//Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound=right_bound_();//Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges=nb_ranges_();//contains the number of ranges
    Vector<Size>& left_bound_index=left_bound_index_();//Specifies the corresponding freq index to the left bounds of the ranges
    Vector<Size>& right_bound_index=right_bound_index_();//Specifies the corresponding freq index to the right bounds of the ranges

    nb_ranges=0;
    //First, we specify which quadratues are counted (every one except the outmost 2 ones depending on discretization direction)
    Vector<Size>& quad_range_weight=quad_range_weight_();
    Size& tot_quad_range_weight=tot_quad_range_weight_();
    tot_quad_range_weight=model.parameters->nquads()-2;//Number of used quadratures is the same for both cases, so we do not need to worry about this
    if (is_upward_disc)
    {//then everything except the final two quads may be used for determining the range
        Size quadidx=0;
        while (quadidx<model.parameters()->nquads()-2)
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }

        while (quadidx<model.parameters()->nquads())
        {
            quad_range_weight[quadidx]=0;
            quadidx++;
        }
    }
    else
    {//then everything except the first two quads may be used for determining the range
        Size quadidx=0;
        while (quadidx<2)
        {
            quad_range_weight[quadidx]=0;
            quadidx++;
        }

        while (quadidx<model.parameters()->nquads())
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }
    }

    //Dummy initialization, as we do not know a priori what the first counted quadrature will be
    Real leftbound=0.0;
    Real rightbound=0.0;
    Size leftboundidx=0;
    Size rightboudidx=0;

    bool leftbound_specified=false;//Specifies whether the left bound for the current range has already been set

    Vector<Size>& line_count=line_count_();//counts the number of quadratures encountered for each line
    //This can get complicated due to having no requirement to use the same quadrature for all lines...
    //Due to having no requirement that all lines should use the same quadrature, we will just try to compute continguous ranges by iterating over the frequencies itself.
    //Thus we iterate over all freqs, and count the occurences of the lines; if we encounter a line the maximal number of times, set line count to zero and check if continguous region has ended
    //new range gets added if no other line is overlapping.
    //This only works because the frequencies themselves are sorted
    for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
    {
        //get corresponding line of the freq
        const Size lineidx=model.radiation.frequencies.corresponding_line_matrix(curr_point, freqidx);
        //FIXME: REPLACE corresponding_z_for_line WITH FUTURE VERSION WHICH ALSO TAKES INTO ACCOUNT THE POINT INDEX
        const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[freqidx];
        const Real curr_freq=model.radiation.frequencies.nu(curr_point, freqidx)*curr_shift;

        if (!leftbound_specified&&(quad_range_weight[quadidx]))
        {//specify left bound if not yet done so (and the point actually belongs to any of our ranges)
            leftbound_specified=true;
            leftbound=curr_freq;
            leftboundidx=freqidx;
        }
        rightbound=curr_freq;
        rightboundidx=freqidx;
        //increment number of freqs encountered // not incrementing the if boundary quadratures are encountered
        line_count[lineidx]+=quad_range_weight[quadidx];
        //if we have counted all relevant quads belonging to a specific line:
        if (line_count[lineidx]==tot_quad_range_weight)
        {
            line_count[lineidx]=0;
            //Then if no other lines are currently being counted (i.e. entire count is zero),
            // we conclude that the current continguous range has ended
            if (std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
            {
                //add range
                left_bound[nb_ranges]=leftbound;
                right_bound[nb_ranges]=rightbound;
                left_bound_index[nb_ranges]=leftboundidx;
                right_bound_index[nb_ranges]=rightboundidx;
                nb_ranges++;

                //set flag for new left bound
                leftbound_specified=false;
            }
        }
    }
    //TODO: check line_count to be uniformly 0!!
}


// //TODO: add the freq derivative bdy conditions here TOO
// //maybe TODO: add some array of bools for determining this? (or just modify the if-clause, checking whether and flag is set TODO ADD AS ARGUMENT and index <)
// //We need the curr_point, as the ray_direction is not given
// //Sets the left boundary frequencies needed due to shift
// inline void Solver :: set_left_implicit_boundary_frequencies(const Size nextpoint, const Size currpoint, const Size lineindex,
//                           std::vector<std::array<std::deque<std::tuple<Real, Size, Size>>>>& tuple_deque_array_vector, bool disc_boundary_condition)
// {//, const Size splitindex
//     //freq splits no longer used; just use all individual lines after checking for overlap?
//     //iterate over the sorted lines?
//     //or use the computed range of the current point?
//     // const Size min_freq_idx=freqsplits_()[splitindex];//freq index for the current point
//     // const Size max_freq_idx=freqsplits_()[splitindex+1]-1;//minimal freq index for this split
//     const Size min_freq_idx=freqsplits_()[splitindex];//freq index for the current point
//     const Size max_freq_idx=freqsplits_()[splitindex+1]-1;//minimal freq index for this split
//     Real curr_min_bound=model.radiation.frequencies.nu(currpoint, min_freq_idx);
//     //? replace by line_quad_right_overlap_(eval)?
//     if (disc_boundary_condition)//just add the outermost two freqs without checking if we decide to have a boundary condition anyway
//     {
//         //FIXME: ordering can be completely wrong? //Yes, this is due to arbitrarily deciding that the outermost two freqs must be a boundary condition
//         //Er, no, the boundary condition issues are already solved in the same iteration, as the innermost freqs will be filled in by the boundary conditions itself (and thus popped from the deque)
//         tuple_deque_array_vector[splitindex][0].push_back(std::make_tuple(freq, nextpoint, min_freq_idx));
//         tuple_deque_array_vector[splitindex][0].push_back(std::make_tuple(freq, nextpoint, min_freq_idx+1));
//         min_freq_idx+=2;//first two indices have now already been processed
//     }
//     //Wait a minute, shouldn't we compare versus the entire previous range?
//     for (Size next_freq_idx=min_freq_idx; next_freq_idx<=max_freq_idx; next_freq_idx++)
//     {
//         const Real temp_next_freq = model.radiation.frequencies.nu(nextpoint, next_freq_idx);
//         if (temp_next_freq>=curr_min_bound)//boundary condition not needed if falls inside range of currpoint
//         {
//             break;
//         }
//         //add at the right of the current boundary frequencies needed
//         tuple_deque_array_vector[splitindex][0].push_back(std::make_tuple(freq, nextpoint, next_freq_idx));
//     }
// }
//
// //TODO: add the freq derivative bdy conditions here TOO
// //maybe TODO: add some array of bools for determining this? (or just modify the if-clause, checking whether and flag is set TODO ADD AS ARGUMENT and index <)
// //We need the curr_point, as the ray_direction is not given
// //Sets the left boundary frequencies needed due to shift
// inline void Solver :: set_left_implicit_boundary_frequencies(const Size nextpoint, const Size currpoint,
//                           std::vector<std::array<std::deque<std::tuple<Real, Size, Size>>>>& tuple_deque_array_vector, bool disc_boundary_condition)
// {
//     //Assumes we have already computed the ranges belonging to currpoint
//     //Also assumes that the overlap between lines of nextpoint has been computed
//     //First: for all lines, set the outermost two points (if non-overlapping and correct direction) as boundary
//
//     //Second: just iterate over all lines (from left to right?), checking whether they lie within the specified ranges (O(nfreqs))
//     //if not, then also add them to their respective
//     TODO COMMENTS
//     //mhh, cant we just iterate over all freqs (is simpler than iterating over all lines, then figuring out what indices overlap)
//     for (ALL FREQS)
//     {
//         //wait a minute, cant we also compute whether overlap exists at the same time?
//         //more than a single line should then be currently counted
//         //disc_boundary_condition=!is_upward_disc
//         if (disc_boundary_condition&&INDEX=0,1&&!overlap)
//         {
//             //certainly boundary condition
//         }
//         //otherwise check versus computed ranges
//     }
// }
//
// public void check_for_boundary()
// {
//
// }

//Sets all implicit boundary conditions, using the computed frequency ranges for currpoint and overlaps for nextpoint
inline void Solver :: set_implicit_boundary_frequencies(const Size nextpoint, const Size currpoint, Real shift_next,
                std::multi_map<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index, bool is_upward_disc)
{
    // //OR just save in vector which lines are overlapping and do (not?) strictly need boundary conditions at nextpoint
    // if (is_upward_disc)
    // {
    //     //TODO: check whether to use sorted lines (should not matter too much due to insert)
    //     for (Size lineidx=0; lineidx<model.parameters->nlines(); lineidx++)
    //     {
    //         if (!overlappingTODO)
    //         {
    //             Real bdyfreq=model.radiation.frequencies.nu(nextpoint, next_freq_idx)*shift_next;
    //             multimap_freq_to_bdy_index.insert(freq, std::make_tuple(nextpoint, next_freq_idx?));
    //         }
    //     }
    // }
    // else
    // {
    //
    // }

    Size curr_range_index=0;
    Vector<Size>& left_bound=left_bound_();
    Vector<Size>& right_bound=right_bound_();
    Vector<unsigned char>& line_quad_overlap=line_quad_discdir_overlap_();
    Size nb_ranges=nb_ranges_();
    Real left_bound_freq=left_bound[curr_range_index];
    Real right_bound_freq=right_bound[curr_range_index];
    //Or just do all at once, checking the properties of each freq sequentially
    //duplicated due to difference due to discretization direction
    if (is_upward_disc)
    {
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            Real next_freq=model.radiation.frequencies.nu(nextpoint, freqidx)*shift_next;
            const Size lineidx=model.radiation.frequencies.corresponding_line[freqidx];
            const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[freqidx];
            // (quadidx==model.parameters->nquads()-1||quadidx==model.parameters->nquads()-2)
            // replaced by quadidx>=nquads()-2
            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx>=nquads()-2&&!line_quad_overlap[lineidx])
            {
                //then is required boundary condition
                multimap_freq_to_bdy_index.insert(next_freq, std::make_tuple(nextpoint, freqidx));
                continue;
            }
            //ranges are sorted, just like the freqs; therefore we can relatively efficently compare the non-overlap
            // by just advancing the index for the ranges
            while(right_bound_freq<next_freq && curr_range_index<nb_ranges-1)
            {
                curr_range_index++;
                Real left_bound_freq=left_bound[curr_range_index];
                Real right_bound_freq=right_bound[curr_range_index];
            }
            //freq outside range if <left_bound || >right_bound
            if (next_freq<left_bound_freq || next_freq>right_bound_freq)
            {
                //then is required boundary condition
                multimap_freq_to_bdy_index.insert(next_freq, std::make_tuple(nextpoint, freqidx));
            }
        }
    }
    else
    {
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            Real next_freq=model.radiation.frequencies.nu(nextpoint, freqidx)*shift_next;
            const Size lineidx=model.radiation.frequencies.corresponding_line[freqidx];
            const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[freqidx];
            // (quadidx==0||quadidx==1)
            // replaced by quadidx<2
            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx<2&&!line_quad_overlap[lineidx])
            {
                //then is required boundary condition
                multimap_freq_to_bdy_index.insert(next_freq, std::make_tuple(nextpoint, freqidx));
                continue;
            }
            //ranges are sorted, just like the freqs; therefore we can relatively efficently compare the non-overlap
            while(right_bound_freq<next_freq && curr_range_index<nb_ranges-1)
            {
                curr_range_index++;
                Real left_bound_freq=left_bound[curr_range_index];
                Real right_bound_freq=right_bound[curr_range_index];
            }
            //freq outside range if <left_bound || >right_bound
            if (next_freq<left_bound_freq || next_freq>right_bound_freq)
            {
                //then is required boundary condition
                multimap_freq_to_bdy_index.insert(next_freq, std::make_tuple(nextpoint, freqidx));
            }
        }
    }
}

// Matches the overlapping boundary conditions, using the currently computed ranges?
inline void Solver :: match_overlapping_boundary_conditions(Model& model, const Size currpoint, const Real curr_shift, std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index, TODO?)
{
    //just iterate simultaneously over the ranges and frequencies
    //TPDO first define the indices corresponding to the ranges//DONE
    Vector<Real>& left_bound=left_bound_();//Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound=right_bound_();//Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges=nb_ranges_();//contains the number of ranges
    Vector<Size>& left_bound_index=left_bound_index_();//Specifies the corresponding freq index to the left bounds of the ranges
    Vector<Size>& right_bound_index=right_bound_index_();//Specifies the corresponding freq index to the right bounds of the ranges

    Size curr_range_index=0;
    //Assumption: at least a single line exists -> a single range exists
    Real curr_range_min=left_bound[curr_range_index];
    Real curr_range_max=right_bound[curr_range_index];
    Real curr_range_freq_idx=curr_range_min;
    //TODO: also add frequencies!
    //sorted iterating over all boundary frequencies
    //TODO: just get iterator over multimap and simulataneously check for end of ranges too!
    std::multimap<Real, std::tuple<Size, Size>>::iterator it=multimap_freq_to_bdy_indexs.begin();
    while (it!=multimap_freq_to_bdy_index.end()&&curr_range_index<nb_ranges)
    {
        const Real bdy_freq=it->first;//get boundary freq
        //compare the frequency to the max freq of current range; if larger, then we need to compare against the next range
        while (right_bound[curr_range_index]<bdy_freq)
        {
            if (curr_range_index==nb_ranges)
            {
                //if curr_range_index==nb_ranges, then all remaining boundary frequencies lie outside of the ranges
                //so we do not need to waste any more time on trying to match them
                return;//
            }
            curr_range_index++;
            //update info about ranges
            Real curr_range_min=left_bound[curr_range_index];
            Real curr_range_max=right_bound[curr_range_index];
            Real curr_range_freq_idx=curr_range_min;
        }
        //now it is guaranteed that rightbound>=bdy_freq, so we only need to check whether leftbound<=bdy_freq
        if (TODO INSIDE RANGE)
        {
            //Do the boundary stuff
            //but first check which freqs are exactly in front or behind the

            //while(&& <) looping
            while (TODO bdy_freq>TODORANGEFREQ&&curr_range_freq_idx<curr_range_max)
            {
                curr_range_freq_idx++;
            }
            //the curr range freq should now be larger/equal to the bdy freq
            //EXCEPT in the case that bdy_freq==freq[curr_range_freq_idx]

            //Now we are finally ready to apply the boundary conditions
            TODO: dIdnu stuff
            //set explicit coefficents to TWICE the normal 1/Δν value

            //set dtau to approx 0

            //set S? (nah just ignore the existence, as dtau≃0)

            //pop value from iterator
            it=multimap_freq_to_bdy_index.erase(it);
        }
        else
        {
            it++;
        }
    }

    for (auto tuple: multimap_freq_to_bdy_index)
    {
        //if the boundary conditions lies within the current range
    }

}


// //TODO: add the freq derivative bdy conditions here TOO
// //We need the curr_point, as the ray_direction is not given
// //Sets the left boundary frequencies needed due to shift
// inline void Solver :: set_right_implicit_boundary_frequencies(const Size nextpoint, const Size currpoint, const Size splitindex,
//                           std::vector<std::array<std::deque<std::tuple<Real, Size, Size>>>>& tuple_deque_array_vector)
// {
//     //Using Size, so add +1 to all indices not do 0-1
//     const Size min_freq_idx=freqsplits_()[splitindex]+1;//freq index for the current point
//     const Size max_freq_idx=freqsplits_()[splitindex+1]-1+1;//minimal freq index for this split
//     Real curr_max_bound=model.radiation.frequencies.nu(currpoint, max_freq_idx-1);
//     if (disc_boundary_condition)//just add the outermost two freqs without checking if we decide to have a boundary condition anyway
//     {
//         //FIXME: ordering can be completely wrong? //No: see comment in function above
//         tuple_deque_array_vector[splitindex][1].push_back(std::make_tuple(freq, nextpoint, max_freq_idx-1));
//         tuple_deque_array_vector[splitindex][1].push_back(std::make_tuple(freq, nextpoint, max_freq_idx-2));
//         max_freq_idx-=2;//last two indices have now already been processed
//     }
//     for (Size next_freq_idx=max_freq_idx; next_freq_idx>=max_freq_idx; next_freq_idx--)
//     {
//         const Real temp_next_freq = model.radiation.frequencies.nu(nextpoint, next_freq_idx-1);
//         if (temp_next_freq<=curr_max_bound)//boundary condition not needed if falls inside range of currpoint
//         {
//             break;
//         }
//         //add at the right of the current boundary frequencies needed
//         tuple_deque_array_vector[splitindex][1].push_front(std::make_tuple(freq, nextpoint, next_freq_idx-1));
//     }
// }



// Maybe TODO? add some conversion from (point,freq) to general index? Only useful if we were not using Matrix
// Mhh, vector of size 2 is sufficient for our pu

//Idea, maybe specialize this function by setting the discretization direction in a template? So less if-clauses should appear
//BEFORE EVERY RAY IN EVERY DIRECTION (also backwards, even though it is seemingly mirrored)
//Run twice; once for forward ray, once for backward ray TODO: implement backward ray version
//assumes ray has already been traced
//?UNLESS STATED OTHERWISE, WE WILL USE THE NON-DOPPLER SHIFTED FREQUENCIES FOR DETERMINING BOUNDARY CONDITIONS
//wait, that does not make sense, we must doppler shift them obviously
inline void Solver :: comoving_ray_bdy_setup(TODO)
{
    // Size n_freqsplits=n_freqsplits_();
    //Trace the ray as usual, not needing anything except the frequencies
    //Start with imposing some initial boundary conditions
    //Note: this assumes forward ray, not backward ray TODO IMPLEMENT OTHER DIRECTION
    Matrix<Real> intensities=intensities_();
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        intensities(first_index, freqid)=boundary_intensity(model, o, model.radiation.frequencies.nu(o, freqid);
    }

    Size first_index=first_();//?
    rayposidx=last_();//ray position index -> point index through nr_()[rayposidx]
    //from first_+1 to last_ index, trace the ray in
    // for (Size rayposidx=first_(); rayposidx<=last_(); pointidx++)//for loop does at least one iteration, which results in oob access if ray goes trough one point
    // std::vector<std::deque<std::map<Real, std::pair<Size, Size>>>> deque_vector;//Holds temporary values for determining which frequency boundary values are needed where
    // deque_vector.resize(TODO_N_SPLITS);//Wait, due to line width increase, technically up to 2*Nsplit deques are needed for storing the boundary points
    //Therefore a fixed size vector is not appropriate
    //However, as we will only traverse the deques sequentially (for performance savings), I see no reason not to add another deque instead of a vector
    // std::deque<std::deque<std::map<Real, std::pair<Size, Size>>>> deque_deque;//Holds temporary values for determining which frequency boundary values are needed where
    // std::deque<std::deque<std::tuple<Real, Size, Size>>> deque_deque;//Holds temporary values for determining which frequency boundary values are needed where
    // std::vector<std::array<std::deque<std::tuple<Real, Size, Size>>>> tuple_deque_array_vector;//Holds temporary values for determining which frequency boundary values are needed where
    // tuple_deque_array_vector.resize(n_freqsplits);

    //Due to complicated handling of boundary conditions, we will for now forgo the use of these hard to use deques, instead opting for a sane map
    //Will up computational demands to O(Nfreqs*ln(Nfreqs)*Npoints), but that can't be helped //note that this is a very pessimistic estimate, as Nfreqs should be replaced with the average number of boundary conditions created
    // for (Size tempidx=0; tempidx<n_freqsplits; tempidx++)
    // {
    //     tuple_deque_array_vector[tempidx].resize(2);//A left and a right
    // }
    //for every freq split, a deque holding maps of the frequencies to the pair of index values (point, freqidx)

    //dummy initialization; just need some references to Reals for get_eta_and_chi
    //Previously defined vectors are not of the right size for me (O(npointsonray) instead of O(nfreqs))
    //So I just compute per point individually (however the compiler should notice and loop unrolling could still happen)
    Real eta_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real chi_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real eta_curr=0.0;
    Real chi_curr=0.0;

    //helper values for computing the freq derivative terms
    Real dfreqsmall=0.0;
    Real dfreqlarge=0.0;

    //POSSIBLE OPTIMIZATION: only query them once, moving ..._curr to ..._next after each rayposidx change
    //But then this should change to a vector holding these values

    std::multimap<Real, std::tuple<Size, Size>> multimap_freq_to_bdy_index;//boundary conditions can overlap, so multimap it is
    //Yes, for determining boundary conditions, we are going backwards; this seems strange, but is a result of easier treatment of the spatial boundary conditions at the end
    while (rayposidx>first_())
    {
        const Size nextpoint=nr_()[rayposidx];
        const Size currpoint=nr_()[rayposidx-1];

        const Real next_shift=TODO;
        const Real curr_shift=TODO;

        //TODO: for less mistakes, redefine rayposidx, rayposidx-1 using different names

        //For the sake of accurately putting boundary conditions, this should point in the correct direction; otherwise we might accidentally put boundary conditions somewhat closer to the middle of a line
        const bool is_upward_disc=(TODO projectedvnext>=TODO projectedvcurr);
        //TODO MAKE/(find in geometry?) (trivial) FUNCTION FOR THIS

        //Compute whether upward disc is needed (just compute )
        //\hat{n}⋅\bar{v}; TODO: check usual convention in sign of doppler shift again

        //somewhere, we should also compute all of S, Δτ, dI/dν stuff
        //this first part might be appropriate, as we can as well overwrite all this stuff later for the boundary conditions
        //err, Δτ is not available until we know which freq matches which other freq
        match_frequency_indices(model, nextpoint, currpoint, rayposidx, rayposidx-1, TODO is_upward_disc)//O(nfreqs)

        //Now compute all default stuff, computing bogus for the non-matched indices (but as these correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Real nextfreq=model.radiation.frequencies.nu(nextpoint, next_freq_idx);
            const Real currfreq=model.radiation.frequencies.nu(currpoint, start_indices_()(rayposidx, next_freq_idx)[1]);
            //technically, we also need to read the point index from this start_indices_ (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to currpoint at this moment in time
            //TODO GET LINE FROM FREQ ID for ll
            get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta_next, chi_next);
            //Possible optimization, only query eta/chi_curr, as the next value is the same either way
            get_eta_and_chi<TODO APPROX>(model, currpoint, ll?, currfreq, eta_curr, chi_curr);
            const Real Snext=eta_next/chi_next;
            const Real Scurr=eta_curr/chi_curr;
            S_next_()(rayposidx, next_freq_idx)=Snext;
            S_curr_()(rayposidx, next_freq_idx)=Scurr;
            // const Real dtau = trap (chi_c[f], chi_n[f], dZ);
            const Real dtau = trap (chi_curr, chi_next, dZ?);
            delta_tau_()(rayposidx, next_freq_idx)=dtau;
            // Compute_S // done
            // Compute_dtau // done
        }
            // compute_freq_der_coeffs //and indices
            // compute for rayposidx and rayposidx-1 !!

            //TODO refactor, bring if clause outside of this loop (adding another for loop if necessary)
        //This part precomputes the frequency derivative coefficients (and corresponding indices), conveniently ignoring the fact that for some outermost line quadrature frequencies, we must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        //The discretization is evidently depends on freq disc direction
        if (is_upward_disc)
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs()-2; next_freq_idx++)
            {
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices
                Size curr_freq_idxp1=(curr_freq_idx+1)%model.parameters->nfreqs();//index+1
                Size curr_freq_idxp2=(curr_freq_idx+2)%model.parameters->nfreqs();//index+2
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                dfreqsmall=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
                            *shift_curr TODO;
                dfreqlarge=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
                            *shift_curr TODO;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=corresp_curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx+1))
                            *shift_next TODO;
                dfreqlarge=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx))
                            *shift_next TODO;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx+2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx+1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
            //and also just set some irrelevant values for the final boundary points, as these should get overwritten?
            //Wait do we then even need to set it here?
            //TODO: figure out whether I need to set some random values here!
        }
        else
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
            {
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices

                //MODULUS FOR 'negative' (more like overflowing) unsigned integer is in general not that easy to compute
                //However, as i know it only goes ever so slightly negative ('-1,-2'): (mod+i)%mod is sufficient for the computation (mod>=1,2)
                Size curr_freq_idxm1=(curr_freq_idx-1+model.parameters->nfreqs())%model.parameters->nfreqs();//index-1
                Size curr_freq_idxm2=(curr_freq_idx-2+model.parameters->nfreqs())%model.parameters->nfreqs();//index-2
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.nu(currpoint, curr_freq_idxm1)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr TODO;
                dfreqlarge=-(model.radiation.frequencies.nu(currpoint, curr_freq_idx)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr TODO;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=corresp_curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next TODO;
                dfreqlarge=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next TODO;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
            //and also just set some irrelevant values for the final boundary points, as these should get overwritten?
            //Wait do we then even need to set it here?
            //TODO: figure out whether I need to set some random values here!
        }

        //Note to self: we can just use all default freq derivative stuff, except for the 2 outermost points on the grid.
        // Because I define every computation in terms of the next point index, I will have no problems with non-existant/too far frequencies, as I have made the computed curr regions slightly smaller
        //Note: it might seem we will try to use oob values, but these should be overwritten anyway
        //CHECK THIS FOR A SIMPLE EXAMPLE!!!

        //now start classifying the boundary points
        // Size bdy_deque_split_index=0;
        // Size next_split_index=0;
        //For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc, curr_shift);//for the current point, we definitely need the exact ranges (fiddled with to ensure enough bdy conditions)
        // and compute the overlap between lines on next_point
        get_overlapping_lines(model, nextpoint, is_upward_disc);//technically, we should compute which lines overlap only in this part
        //?WHILE LOOP SEEMS WEIRD?

        //First, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be treated as boundary points
        set_implicit_boundary_frequencies(nextpoint, currpoint, shift_next, multimap_freq_to_bdy_index, is_upward_disc);
        // And now we check what boundary conditions need to be evaluated using the intensities at curr_point


        for (ALL LINES)//? my code already runs over all frequencies?
        // for (ALL SPLITS)
        while (bdy_deque_split_index<n_freqsplits &&next_split_index<n_freqsplits)
        {
            //Wait a bit, this left deque can only be fed by going to the right, so we do not need these silly minima for adding stuff to the deque
            // Real left_deque_min=std::get<0>(deque_vector[bdy_deque_split_index][0].front());
            // Real left_deque_max=std::get<0>(deque_vector[bdy_deque_split_index][0].back());
            //While the right deque can only be fed by the outermost freq going to the left!
            //adds the oob freqs to the deque[splitindex]
            set_left_implicit_boundary_frequencies(TODO, using computed ranges)//
            set_right_implicit_boundary_frequencies(TODO, using computed ranges)


            // NOTE TO SELF: already included in boundary conditions
            // //also set the 2 last ones in the discretization direction to be boundary points //probably far enough away from line center
            // //TODO CHECK REASONABLE MIN freq dist; (optical depth increments will also limit this approach, as in case of near infinite optical depths, the affected frequency range can be arbitrarily large)
            // if (forwardfreqdisc)
            // {
            //     std::deque<std::tuple<Real, Size, Size>>& left_deque=tuple_deque_array_vector[bdy_deque_split_index][0];
            //     //check whether the two rightmost things are included
            //     //just grab two last indices if possible
            //     if (left_deque.size()>=2)
            //     {
            //
            //     }
            //     else if (left_deque.size()==1)
            //     {
            //
            //     }
            // }
            // else
            // {
            //
            // }
            // //add
        }

        bdy_deque_split_index=0;//TODO replace by line index, as we will interpolate using line quadratures. (makes me not have to deal with defining some suitable distance measure between lines)
        Size curr_split_index=0;
        while (bdy_deque_split_index<n_freqsplits && curr_split_index<n_freqsplits)
        {//mashing two loops together, saving one order of magnitude of time
        // for (ALL SPLITS (saved boundary))//iterating over all boundary splits (sorted obviously)
        // {//TODO TODO MASH THESE LOOPS TOGETHER, saving one order of magnitude (N_lines) of time
        //     for (ALL SPLITS (current splits))//iterating over all current splits (also sorted)
            //States that the matching (overlapping with the explicit part) boundary frequencies, should be removed from the dequeus, filling in the boundary conditions
            //Evidently iterate over all non-empty deques
            //NOTE TO SELF: for the boundary conditions, set Δτ very low, the implicit derivative coeffs to 0, the explicit derivative coeffs to TWICE the usual coeffs (as normally a part of the freq derivative is implicit).
            //and the source function can be ignored (Δτ very low)
            match_implicit_boundary_frequencies
            //only at the very end (initial boundary conditions), the source must be set very high and the optical depth too. (then also set freq der coeffs to zero)
        // }
        }
        //if tracing the ray forward for discovering the boundary conditions
        //get curr point oob freqs
        //store ?map? in deque (also store limits (or just query them)) (map should map freq to exact index (pair?) (point, freqidx) which put it on)

        //get next point oob freqs

        //OR if tracing the ray backward to discover the boundary conditions
        //then store 'next'? point oob freqs in deque and also their indices (so map again necessary// err, is exactly the same procedure, so just choose one

        //AT THE END, overwrite the lastmost two freq derivative
    }



    //For every point on the ray
    for(TODO)
    {
        //'save' curr oob intensity, overwriting if necessary (easy way: if previous iteration stored and too close (0.? line widths))
        //err, overwriting is a bit unnecessary in this step, as it should only be deleted when consumed to determine an intensity
        //? should we also consider the furthermost intensity as a boundary condition ?
        //
        // So the

        //for the next point, check if oob intensity lies close enough between 2 values in the list
        //Also include curr point farmost freq in search?
        //By default, bdy condition; but if one nearby stored frequency is found, one can also check the split for previous freqs to fully enclose it
        //if not found, evidently we need initial intensities
    }
    //Also compute the freq derivative coeffs here? Depends on the discretization direction, so does not strictly fit in the same piece of code
    //However, this would lead us not to duplicate logic..., thus might be useful
    //If we use this opportunity, we can save some headaches



    //In the same vein, cant we also precompute all Δτ, S? (same memory requirement)
    //Then the boundary condition can be obtained by cheating (set Δτ=very high (such that e^-Δτ≃0), S=bdy intensity), intuitively this is just rethermalizing the intensity
    //This is a bit optimistic, as we need S for two points... (however duplicating S is possible (Scurr and Snext for every position increment)
    //; so we might be able to choose the intensity by only changing a single increment
    //
    //Then the only thing left would be the (quite simplified) computation

    //And even save the traced points, as we do not really need to trace the same ray once again...
    //However, can't we jsut give the ray points to this function? (just vector with points) //should also be able to be reversed without performance penalty

}




//
inline void Solver :: compute_freq_der_term(TODO)
{
    //Assumes we have already computed the giant lists of boundary freqs?

}



// TODO: CHANGE THIS TO COMOVING SOLVER
// Note: single line approximation might be useful for determining the splits. IMPLEMENT FANCY VERSION
template<ApproximationType approx>
inline void Solver :: solve_comoving_order_2_sparse (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        accelerated_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }

    // For each ray, solve the radiative transfer equation in a comoving manner
    accelerated_for (rr, model.parameters->hnrays(),
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        // const Real     wt = model.geometry.rays.weight   [rr] * two;
        // const Vector3D nn = model.geometry.rays.direction[rr];

        std::cout << "--- rr = " << rr << std::endl;

        //for every ray to trace
        const Size n_rays_to_trace=points_to_trace_ray_through.size();
        accelerated_for (rayidx, n_rays_to_trace,
        {
            //trace and solve over the ray in both directions, incrementing J as necessary
            //too complicated to decouple ray tracing from solving
            solve_comoving_order_2_sparse(model, o, rr, rayidx, dshift_max);
            // solve_comoving_order_2_sparse(model, o, ar, rayidx, dshift_max); //reuse of some data is possible
        })
    })

    // TODO REPLACE/DELETE THIS?
    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr] * two;
        const Vector3D nn = model.geometry.rays.direction[rr];

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);

                nr_   ()[centre] = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                n_tot_() = (last_()+1) - first_();

                if (n_tot_() > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o,  lspec.nr_line[o][k][z]);

                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda           (model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}

//
inline void Solver :: solve_comoving_order_2_sparse (
      Model& model,
      const Size o,//ray origin point
      const Size r,//ray direction index
      const Size rayidx//ray trace index
      const double dshift_max)
{
    Vector<Real>& eta_c = eta_c_();//current emissivity for all freqs
    Vector<Real>& eta_n = eta_n_();//next emissivity    for all freqs

    Vector<Real>& chi_c = chi_c_();//current opacity    for all freqs
    Vector<Real>& chi_n = chi_n_();//next opacity       for all freqs

    Vector<Real>& curr_intensity = curr_intensity_();//current intensity for all freqs
    Vector<Real>& next_intensity = next_intensity_();//storage for next intensity for all freqs

    Vector<unsigned char>& real_pt= real_pt_();//for denoting whether the given point is a real point

    Vector<double>& shift = shift_();// contains the doppler shifts for all points on the ray

    const Real     wt = model.geometry.rays.weight   [rr];//ray direction weight

    double Z = 0.0; //distance along ray
    double dZ= 0.0; //last distance increment

    //trace the ray
    first_() = trace_ray_indicate_point <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    last_ () = trace_ray_indicate_point <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;

    nr_   ()[centre] = o;
    shift_()[centre] = 1.0;
    n_tot_() = (last_()+1) - first_();

    //Now is the perfect time to setup the boundary conditions (and dIdnu stuff?)


    //assuming boundary points are real points ofcourse

    //Set the boundary conditions, as this is the first point on the ray
    //TODO PROBABLY NO LONGER USEFUL HERE
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        curr_intensity[freqid]=TODOBDY
    }

    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[first_()]))
    {
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];
            lspec.J(nr[first_()],k) += lspec.quadrature.weights[z] * wt * curr_intensity[freqid];// Su_()[centre];
        }
    }

    rayposidx=first_()+1;//ray position index -> point index through nr[rayposidx]
    //from first_+1 to last_ index, trace the ray in
    // for (Size rayposidx=first_(); rayposidx<=last_(); pointidx++)//for loop does at least one iteration, which results in oob access if ray goes trough one point
    while (rayposidx<=last_())
    {
        //Do the freq split
        const Size currpointidx=nr[rayposidx-1];
        const Size nextpointidx=nr[rayposidx];

        split_frequencies(model, nextpointidx);

        const Real rel_doppler_shift=shift[rayposidx]-shift[rayposidx-1];

        if (rel_doppler_shift>0)
        {
            //For every split, do freq mismatch and get boundary conditions
            //also rotate the frequencies if necessary before overwriting with boundary conditions
            for (Size splitid=0; split<n_freqsplits_(); split++)
            {
                //determine number of boundary conditions
                const Size nboundaryconditions=std::max(2,TODO);//at least two are needed for the solver

                // RESET frequency correspondence?
                //For all other (next) frequencies, determine what previous frequency corresponds
                // and 'shift' the current intensity accordingly



              //set boundary stuff depending on sign doppler shift


              //get explicit term frequency derivative term (for all but the outermost two points)
              //and do the explicit part


            }

            //Do freq mismatch and get boundary conditions?

            //for every split
            //Actually compute the intensity (explicit part)

            //for every split
            //Actually compute the intensity (implicit part)

            //set Boundary conditions? (or do this when getting bdy conditions)

        }
        else
        {//almost exactly the same, but in reverse direction
            //TODO
        }


        //Finally increment J if and real point check if closest ray
        if (real_pt[rayposidx]&&closest_ray(rr, nextpointidx))
        {
            //then obviously add (weighted) to J
            for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
            {
                // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
                // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
                // Size1 corresponding_z_for_line;           ///< quadrature number corresponding to frequency
                const Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];
                const Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];
                lspec.J(nextpointidx,k) += lspec.quadrature.weights[z] * wt * curr_intensity[freqid];// Su_()[centre];
            }
        }

        //and finally move the contents of next_intensity to curr_intensity
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            curr_intensity[freqid]=next_intensity[freqid];
        }

        rayposidx++;
    }





}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_anis (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.J      .resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_0   .resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J       (o,k) = 0.0;
                lspec.J2_0    (o,k) = 0.0;
                lspec.J2_1_Re (o,k) = 0.0;
                lspec.J2_1_Im (o,k) = 0.0;
                lspec.J2_2_Re (o,k) = 0.0;
                lspec.J2_2_Im (o,k) = 0.0;
            }
        })
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr];
        const Vector3D nn = model.geometry.rays.direction[rr];

        const Real wt_0    =     inv_sqrt2 * (three * nn.z() * nn.z() - one);
        const Real wt_1_Re =        -sqrt3 *          nn.x() * nn.z();
        const Real wt_1_Im =        -sqrt3 *          nn.y() * nn.z();
        const Real wt_2_Re =  half * sqrt3 * (nn.x() * nn.x() - nn.y() * nn.y());
        const Real wt_2_Im =         sqrt3 *          nn.x() * nn.y();

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);

                nr_   ()[centre] = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                n_tot_() = (last_()+1) - first_();

                if (n_tot_() > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o, lspec.nr_line[o][k][z]);

                            const Real du = lspec.quadrature.weights[z] * wt * Su_()[centre];

                            lspec.J       (o,k) += two     * du;
                            lspec.J2_0    (o,k) += wt_0    * du;
                            lspec.J2_1_Re (o,k) += wt_1_Re * du;
                            lspec.J2_1_Im (o,k) += wt_1_Im * du;
                            lspec.J2_2_Re (o,k) += wt_2_Re * du;
                            lspec.J2_2_Im (o,k) += wt_2_Im * du;
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            const Real du = lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));

                            lspec.J       (o,k) += two     * du;
                            lspec.J2_0    (o,k) += wt_0    * du;
                            lspec.J2_1_Re (o,k) += wt_1_Re * du;
                            lspec.J2_1_Im (o,k) += wt_1_Im * du;
                            lspec.J2_2_Re (o,k) += wt_2_Re * du;
                            lspec.J2_2_Im (o,k) += wt_2_Im * du;
                        }
                    }
                }
            })
        }
    }
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2 (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize (                           model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    distributed_for (rr, rr_loc, model.parameters->hnrays(),
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    solve_feautrier_order_2 <approx> (model, o, f);

                    model.radiation.u(rr_loc,o,f)  = Su_()[centre];
                    model.radiation.J(       o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];

                    update_Lambda (model, rr, f);
                }
            }
            else
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    model.radiation.u(rr_loc,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.J(       o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr_loc,o,f);
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
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.MPI_gather();}
    cout << "Done MPI gathering Lambda." << endl;
}


inline void Solver :: image_feautrier_order_2 (Model& model, const Size rr)
{
    Image image = Image(model.geometry, Intensity, rr);

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters->npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image_feautrier_order_2 (model, o, f);

                image.I(o,f) = two*Su_()[last_()] - boundary_intensity(model, nr_()[last_()], model.radiation.frequencies.nu(o, f));
            }
        }
        else
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image.I(o,f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back (image);
}


inline void Solver :: image_feautrier_order_2_for_point (Model& model, const Size rr, const Size p)
{
    // Redefine p to keep everything as similar to image_feautrier_order_2 as possible
    const Size o = p;


    const Size ar = model.geometry.rays.antipod[rr];

    const Real dshift_max = get_dshift_max (model, o);

    nr_   ()[centre] = o;
    shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

    first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
    n_tot_() = (last_()+1) - first_();

    model. eta_ray.resize (n_tot_(), model.parameters->nfreqs());
    model. chi_ray.resize (n_tot_(), model.parameters->nfreqs());
    model.dtau_ray.resize (n_tot_(), model.parameters->nfreqs());
    model.   u_ray.resize (n_tot_(), model.parameters->nfreqs());

    if (n_tot_() > 1)
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            image_feautrier_order_2_for_point_loc (model, o, f);
        }
    }

}


inline void Solver :: image_optical_depth (Model& model, const Size rr)
{
    Image image = Image(model.geometry, OpticalDepth, rr);

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters->npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image_optical_depth (model, o, f);

                image.I(o,f) = optical_depth_();
            }
        }
        else
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image.I(o,f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back (image);
}


template <Frame frame>
accel inline Size Solver :: trace_ray (
    const Geometry& geometry,
    const Size      o,
    const Size      r,
    const double    dshift_max,
    const int       increment,
          Size      id1,
          Size      id2 )
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z  );

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, nxt, Z, dZ);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z    );

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

//For the comoving solvers, we compute J for all points on a ray,
// so we need to know which points are real
template <Frame frame>
accel inline Size Solver :: trace_ray_indicate_point (
    const Geometry& geometry,
    const Size      o,
    const Size      r,
    const double    dshift_max,
    const int       increment,
          Size      id1,
          Size      id2 )
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    //set self also on ray? No, will be added manually in the middle

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z  );

        set_data_indicate_point (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, crt, Z, dZ);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z    );

            set_data_indicate_point (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}


accel inline void Solver :: set_data (
    const Size   crt,
    const Size   nxt,
    const double shift_crt,
    const double shift_nxt,
    const double dZ_loc,
    const double dshift_max,
    const int    increment,
          Size&  id1,
          Size&  id2 )
{
    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs (dshift);

    if (dshift_abs > dshift_max) // If velocity gradient is not well-sampled enough
    {
        // Interpolate velocity gradient field
        const Size        n_interpl = dshift_abs / dshift_max + 1;
        const Size   half_n_interpl = 0.5 * n_interpl;
        const double     dZ_interpl =     dZ_loc / n_interpl;
        const double dshift_interpl =     dshift / n_interpl;

        if (n_interpl > 10000)
        {
            printf ("ERROR (n_intpl > 10 000) || (dshift_max < 0, probably due to overflow)\n");
        }

        // Assign current cell to first half of interpolation points
        for (Size m = 1; m < half_n_interpl; m++)
        {
            nr   [id1] = crt;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }

        // Assign next cell to second half of interpolation points
        for (Size m = half_n_interpl; m <= n_interpl; m++)
        {
            nr   [id1] = nxt;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }
    }

    else
    {
        nr   [id1] = nxt;
        shift[id1] = shift_nxt;
        dZ   [id2] = dZ_loc;

        id1 += increment;
        id2 += increment;
    }
}

// Comoving solver also needs to know which points are real points in the grid, to be able to compute J from the correct intensity
accel inline void Solver :: set_data_indicate_point (
    const Size   crt,
    const Size   nxt,
    const double shift_crt,
    const double shift_nxt,
    const double dZ_loc,
    const double dshift_max,
    const int    increment,
          Size&  id1,
          Size&  id2 )
{
    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();
    Vector<unsigned char>& real_pt= real_pt_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs (dshift);

    if (dshift_abs > dshift_max) // If velocity gradient is not well-sampled enough
    {
        // Interpolate velocity gradient field
        const Size        n_interpl = dshift_abs / dshift_max + 1;
        const Size   half_n_interpl = 0.5 * n_interpl;
        const double     dZ_interpl =     dZ_loc / n_interpl;
        const double dshift_interpl =     dshift / n_interpl;

        if (n_interpl > 10000)
        {
            printf ("ERROR (n_intpl > 10 000) || (dshift_max < 0, probably due to overflow)\n");
        }

        Size startid1=id1;
        // Assign current cell to first half of interpolation points
        for (Size m = 1; m < half_n_interpl; m++)
        {
            nr   [id1] = crt;
            real_pt[id1]= false;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }

        // Assign next cell to second half of interpolation points
        for (Size m = half_n_interpl; m <= n_interpl; m++)
        {
            nr   [id1] = nxt;
            real_pt[id1]= false;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }
        real_pt[id1-increment]=true;
    }

    else
    {
        nr   [id1] = nxt;
        real_pt[id1]= true;
        shift[id1] = shift_nxt;
        // std::cout<<"shift: "<<shift_nxt<<std::endl;
        dZ   [id2] = dZ_loc;

        id1 += increment;
        id2 += increment;
    }
}


///  Gaussian line profile function
///    @param[in] width : profile width
///    @param[in] diff  : frequency difference with line centre
///    @return profile function evaluated with this frequency difference
////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: gaussian (const Real inverse_width, const Real diff) const
{
    const Real sqrt_exp = inverse_width * diff;

    return inverse_width * INVERSE_SQRT_PI * exp (-sqrt_exp*sqrt_exp);
}


///  Planck function
///    @param[in] temp : temperature of the corresponding black body
///    @param[in] freq : frequency at which to evaluate the function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: planck (const Real temp, const Real freq) const
{
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1 (HH_OVER_KB*freq/temp);
}


///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] p      : point index of the boundary point
///    @param[in] freq   : frequency at which to evaluate boundary condition
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: boundary_intensity (const Model& model, const Size p, const Real freq) const
{
    const Size bdy_id = model.geometry.boundary.point2boundary[p];

    switch (model.geometry.boundary.boundary_condition[bdy_id])
    {
        case Zero    : return 0.0;
        case Thermal : return planck (model.geometry.boundary.boundary_temperature[bdy_id], freq);
        default      : return planck (T_CMB, freq);
    }
}


///  Getter for the emissivity (eta) and the opacity (chi)
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template<>
accel inline void Solver :: get_eta_and_chi <None> (
    const Model& model,
    const Size   p,
    const Size   ll,  // dummy variable
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    // Initialize
    eta = 0.0;
    chi = 1.0e-26;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters->nlines(); l++)
    {
        const Real diff = freq - model.lines.line[l];
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);
    }
}


///  Getter for the emissivity (eta) and the opacity (chi) in the "one line" approixmation
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  l     : line index corresponding to the frequency
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////////////////////////////////////
template<>
accel inline void Solver :: get_eta_and_chi <OneLine> (
    const Model& model,
    const Size   p,
    const Size   l,
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    const Real diff = freq - model.lines.line[l];
    const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

    eta = prof * model.lines.emissivity(p, l);
    chi = prof * model.lines.opacity   (p, l) + 1.0e-26;
}


///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap (const Real x_crt, const Real x_nxt, const double dZ)
{
    return half * (x_crt + x_nxt) * dZ;
}





accel inline void Solver :: solve_shortchar_order_0 (
          Model& model,
    const Size   o,
    const Size   r,
    const double dshift_max)
{
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& tau = tau_();


    double  Z = 0.0;   // distance along ray
    double dZ = 0.0;   // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point (nxt))
    {
        double shift_c = 1.0;
        double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);
            const Size l    = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi <None> (model, crt, l, freq,         eta_c[f], chi_c[f]);
            get_eta_and_chi <None> (model, nxt, l, freq*shift_n, eta_n[f], chi_n[f]);

            const Real drho = trap (eta_c[f], eta_n[f], dZ);
            const Real dtau = trap (chi_c[f], chi_n[f], dZ);

            tau[f]                   = dtau;
            model.radiation.I(r,o,f) = drho * expf(-tau[f]);
        }

        while (model.geometry.not_on_boundary (nxt))
        {
            crt     = nxt;
            shift_c = shift_n;
              eta_c =   eta_n;
              chi_c =   chi_n;

            model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_n);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                const Real freq = model.radiation.frequencies.nu(o, f);
                const Size l    = model.radiation.frequencies.corresponding_line[f];

                get_eta_and_chi <None> (model, nxt, l, freq*shift_n, eta_n[f], chi_n[f]);

                const Real drho = trap (eta_c[f], eta_n[f], dZ);
                const Real dtau = trap (chi_c[f], chi_n[f], dZ);

                tau[f]                   += dtau;
                model.radiation.I(r,o,f) += drho * expf(-tau[f]);
            }
        }

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f) += boundary_intensity(model, nxt, freq*shift_n) * expf(-tau[f]);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
        }
    }

    else
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f)  = boundary_intensity(model, crt, freq);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
        }
    }
}


accel inline void Solver :: update_Lambda (Model &model, const Size rr, const Size f)
{
    const Frequencies    &freqs     = model.radiation.frequencies;
    const Thermodynamics &thermodyn = model.thermodynamics;

    if (freqs.appears_in_line_integral[f])
    {
        const Size first = first_();
        const Size last  = last_ ();
        const Size n_tot = n_tot_();

        Vector<Size  >& nr          = nr_         ();
        Vector<double>& shift       = shift_      ();
        Vector<Real  >& L_diag      = L_diag_     ();
        Matrix<Real  >& L_upper     = L_upper_    ();
        Matrix<Real  >& L_lower     = L_lower_    ();
        Vector<Real  >& inverse_chi = inverse_chi_();

        const Real w_ang = two * model.geometry.rays.weight[rr];

        const Size l = freqs.corresponding_l_for_spec[f];   // index of species
        const Size k = freqs.corresponding_k_for_tran[f];   // index of transition
        const Size z = freqs.corresponding_z_for_line[f];   // index of quadrature point

        LineProducingSpecies &lspec = model.lines.lineProducingSpecies[l];

        const Real freq_line = lspec.linedata.frequency[k];
        const Real invr_mass = lspec.linedata.inverse_mass;
        const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;

        Real frq = freqs.nu(nr[centre], f) * shift[centre];
        Real phi = thermodyn.profile(invr_mass, nr[centre], freq_line, frq);
        Real L   = constante * frq * phi * L_diag[centre] * inverse_chi[centre];

        lspec.lambda.add_element(nr[centre], k, nr[centre], L);

        for (long m = 0; (m < n_off_diag) && (m+1 < n_tot); m++)
        {
            if (centre >= first+m+1) // centre-m-1 >= first
            {
                const long n = centre-m-1;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_lower(m,n) * inverse_chi[n];

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }

            if (centre+m+1 <= last) // centre+m+1 < last
            {
                const long n = centre+m+1;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_upper(m,n) * inverse_chi[n];

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }
        }
    }
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: solve_feautrier_order_2 (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();


    // Get optical properties for first two elements
    get_eta_and_chi <approx> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <approx> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    inverse_chi[first  ] = 1.0 / chi_c;
    inverse_chi[first+1] = 1.0 / chi_n;

    term_c = eta_c * inverse_chi[first  ];
    term_n = eta_n * inverse_chi[first+1];
    dtau_n = half * (chi_c + chi_n) * dZ[first];

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

            C[first] = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = 1.0 / C[first];   // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        // Get new radiative properties
        get_eta_and_chi <approx> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        inverse_chi[n+1] = 1.0 / chi_n;

        term_n = eta_n * inverse_chi[n+1];
        dtau_n = half * (chi_c + chi_n) * dZ[n];

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    if (n_off_diag == 0)
    {
        if (centre < last)
        {
            /// Write economically: G[last] = (B[last] - A[last]) / A[last];
            GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
            GP[last] = GG[last] / (one + GG[last]);

            for (long n = last-1; n > centre; n--) // use long in reverse loops!
            {
                Su[n] += Su[n+1] * FI[n];

                GG[n] = (C[n] * GP[n+1] + one) * inverse_A[n];
                GP[n] = GG[n] / (one + GG[n]);
            }

            Su    [centre] += Su[centre+1] * FI[centre];
            L_diag[centre]  = inverse_C[centre] / (FF[centre] + GP[centre+1]);
        }
        else
        {
            L_diag[centre] = (one + FF[centre-1]) / (Bl_min_Al + Bl*FF[centre-1]);
        }
    }
    else
    {
        /// Write economically: G[last] = (B[last] - A[last]) / A[last];
        GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
        GI[last] = one / (one + GG[last]);
        GP[last] = GG[last] * GI[last];

        L_diag[last] = (one + FF[last-1]) / (Bl_min_Al + Bl*FF[last-1]);

        for (long n = last-1; n > first; n--) // use long in reverse loops!
        {
            Su[n] += Su[n+1] * FI[n];

            GG[n] = (C[n] * GP[n+1] + one) * inverse_A[n];
            GI[n] = one / (one + GG[n]);
            GP[n] = GG[n] * GI[n];

            L_diag[n] = inverse_C[n] / (FF[n] + GP[n+1]);
        }

        Su    [first] += Su[first+1] * FI[first];
        L_diag[first]  = (one + GG[first+1]) / (Bf_min_Cf + Bf*GG[first+1]);

        for (long n = last-1; n >= first; n--) // use long in reverse loops!
        {
            L_upper(0,n+1) = L_diag[n+1] * FI[n  ];
            L_lower(0,n  ) = L_diag[n  ] * GI[n+1];
        }

        for (Size m = 1; (m < n_off_diag) && (m < n_tot-1); m++)
        {
            for (long n = last-1-m; n >= first; n--) // use long in reverse loops!
            {
                L_upper(m,n+m+1) = L_upper(m-1,n+m+1) * FI[n    ];
                L_lower(m,n    ) = L_lower(m-1,n    ) * GI[n+m+1];
            }
        }
    }
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
accel inline void Solver :: image_feautrier_order_2 (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();


    // Get optical properties for first two elements
    get_eta_and_chi <None> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <None> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    inverse_chi[first  ] = 1.0 / chi_c;
    inverse_chi[first+1] = 1.0 / chi_n;

    term_c = eta_c * inverse_chi[first  ];
    term_n = eta_n * inverse_chi[first+1];
    dtau_n = half * (chi_c + chi_n) * dZ[first];

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        // Get new radiative properties
        get_eta_and_chi <None> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        inverse_chi[n+1] = 1.0 / chi_n;

        term_n = eta_n * inverse_chi[n+1];
        dtau_n = half * (chi_c + chi_n) * dZ[n];

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);


    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    // for (long n = last-1; n > first; n--) // use long in reverse loops!
    // {
    //     Su[n] += Su[n+1] * FI[n];
    // }

    // Su[first] += Su[first+1] * FI[first];
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
accel inline void Solver :: image_feautrier_order_2_for_point_loc (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();


    // Get optical properties for first two elements
    get_eta_and_chi <None> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <None> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    inverse_chi[first  ] = 1.0 / chi_c;
    inverse_chi[first+1] = 1.0 / chi_n;

    term_c = eta_c * inverse_chi[first  ];
    term_n = eta_n * inverse_chi[first+1];
    dtau_n = half * (chi_c + chi_n) * dZ[first];

    model. chi_ray(0, f) =  chi_c;
    model. chi_ray(1, f) =  chi_n;

    model. eta_ray(0, f) =  eta_c;
    model. eta_ray(1, f) =  eta_n;

    model.dtau_ray(0, f) = 0.0;
    model.dtau_ray(1, f) = dtau_n;


    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        // Get new radiative properties
        get_eta_and_chi <None> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        inverse_chi[n+1] = 1.0 / chi_n;

        term_n = eta_n * inverse_chi[n+1];
        dtau_n = half * (chi_c + chi_n) * dZ[n];

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        model. chi_ray(n+1-first, f) =  chi_n;
        model. eta_ray(n+1-first, f) =  eta_n;
        model.dtau_ray(n+1-first, f) = dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    model.u_ray(last-first, f) = Su[last];


    for (long n = last-1; n > first; n--) // use long in reverse loops!
    {
        Su[n] += Su[n+1] * FI[n];

        model.u_ray(n-first, f) = Su[n];
    }

    Su[first] += Su[first+1] * FI[first];

    model.u_ray(0, f) = Su[first];
}

/// Image the optical depth
///////////////////////////////////////////////////////////////////////
accel inline void Solver :: image_optical_depth (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Real eta_c, chi_c;
    Real eta_n, chi_n;

    Real tau = 0.0;


    // Get optical properties for first two elements
    get_eta_and_chi <None> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <None> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    tau += half * (chi_c + chi_n) * dZ[first];


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        eta_c = eta_n;
        chi_c = chi_n;

        // Get new radiative properties
        get_eta_and_chi <None> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        tau += half * (chi_c + chi_n) * dZ[n];
    }

    optical_depth_() = tau;
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: solve_feautrier_order_2_uv (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();


    // Get optical properties for first two elements
    get_eta_and_chi <approx> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <approx> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    inverse_chi[first  ] = one / chi_c;
    inverse_chi[first+1] = one / chi_n;

    term_c = eta_c * inverse_chi[first  ];
    term_n = eta_n * inverse_chi[first+1];
    dtau_n = half * (chi_c + chi_n) * dZ[first];

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

            C[first] = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = one / C[first];   // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Sv[first]  = two * inverse_dtau_f * (I_bdy_f - term_c);

    Su[first] /= Bf;
    Sv[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        // Get new radiative properties
        get_eta_and_chi <approx> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        inverse_chi[n+1] = one / chi_n;

        term_n = eta_n * inverse_chi[n+1];
        dtau_n = half * (chi_c + chi_n) * dZ[n];

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
        Sv[n] = (A[n] * Sv[n-1]        ) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Sv[last] = two * inverse_dtau_l * (I_bdy_l - term_n);

    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;
    Sv[last] = (A[last] * Sv[last-1]           ) * (one + FF[last-1]) * denominator;

    if (centre < last)
    {
        for (long n = last-1; n >= centre; n--) // use long in reverse loops!
        {
            Su[n] += Su[n+1] * FI[n];
            Sv[n] += Sv[n+1] * FI[n];
        }
    }

}


accel inline void Solver :: set_eta_and_chi (Model& model, const Size rr) const
{
    model.eta.resize (model.parameters->npoints(), model.parameters->nfreqs());
    model.chi.resize (model.parameters->npoints(), model.parameters->nfreqs());

    for (Size p = 0; p < model.parameters->npoints(); p++)
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            // Extract the Doppler shift
            const double shift = model.geometry.get_shift <Rest> (0, rr, p, 0.0);
            const Real   freq  = model.radiation.frequencies.nu(0, f);
            const Size   l     = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi <None> (model, p, l, freq*shift, model.eta(p,f), model.chi(p,f));
        }
    }
}


accel inline void Solver :: set_boundary_condition (Model& model) const
{
    model.boundary_condition.resize (model.parameters->nboundary(), model.parameters->nfreqs());

    for (Size b = 0; b < model.parameters->nboundary(); b++)
    {
        const Size p = model.geometry.boundary.boundary2point[b];

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(0, f);

            model.boundary_condition(b,f) = boundary_intensity (model, p, freq);
        }
    }
}


inline void Solver :: set_column (Model& model) const
{
    model.column.resize (model.parameters->nrays(), model.parameters->npoints());

    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            model.column(rr, o) = get_column(model, o, rr);
            model.column(ar, o) = get_column(model, o, ar);
        })
    }

}


accel inline Real Solver :: get_column (const Model& model, const Size o, const Size r) const
{
    Real column = 0.0;

    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt))
    {
        Size crt = o;

        column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;

        while (model.geometry.not_on_boundary(nxt))
        {
            crt = nxt;
            nxt = model.geometry.get_next (o, r, nxt, Z, dZ);

            column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;
        }
    }

    return column;
}
