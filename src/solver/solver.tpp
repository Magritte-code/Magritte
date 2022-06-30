#include <map>
#include <tuple>

template <Frame frame>
inline void Solver :: setup (Model& model)
{
    const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size  width = model.parameters->nfreqs();
    const Size  n_o_d = model.parameters->n_off_diag;

    setup (length, width, n_o_d);
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


inline void Solver :: setup (const Size l, const Size w, const Size n_o_d)
{
    length     = l;
    centre     = l/2;
    width      = w;
    n_off_diag = n_o_d;

    //FIXME: check that all frequencies are different!!
    // if (deltafreq==0.0)
    // {   TODO MOVE TO INITIAL GENERAL SETUP (before tracing any rays)
    //     throw runtime_error("Division by zero due to two frequencies being exactly the same");
    // }

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

//l=max length of ray, w=nfreqs, n_o_d is number of off-diagonals; will be ignored
// inline void Solver :: setup_comoving (Model& model, const Size l, const Size w)
inline void Solver :: setup_comoving (Model& model)
{
    length     = 2 * get_ray_lengths_max <Rest> (model) + 1;
    centre     = length/2;
    width      = model.parameters->nfreqs();

    std::cout<<"setup length: "<<length<<std::endl;
    std::cout<<"setup width: "<<width<<std::endl;


    points_to_trace_ray_through.resize(model.parameters->hnrays());
    for (Size i = 0; i < model.parameters->hnrays(); i ++)
    {
        points_to_trace_ray_through[i].resize(model.parameters->npoints());
    }

    //
    n_rays_through_point.resize(model.parameters->hnrays(), model.parameters->npoints());
    min_ray_distsqr     .resize(model.parameters->hnrays(), model.parameters->npoints());
    closest_ray         .resize(model.parameters->hnrays(), model.parameters->npoints());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        // eta_c_       (i).resize (width);
        // eta_n_       (i).resize (width);
        //
        // chi_c_       (i).resize (width);
        // chi_n_       (i).resize (width);
        //
        // tau_         (i).resize (width);

        real_pt_     (i).resize (length);
        start_indices_(i).resize(length, width);
        for (Size j=0; j<length; j++)
        {
            for (Size k=0; k<width; k++)
            {
                start_indices_(i)(j,k).resize(2);
            }
        }

        line_quad_discdir_overlap_(i).resize(model.parameters->nlines());
        left_bound_(i).resize(model.parameters->nlines());
        right_bound_(i).resize(model.parameters->nlines());
        left_bound_index_(i).resize(model.parameters->nlines());
        right_bound_index_(i).resize(model.parameters->nlines());

        line_count_(i).resize(model.parameters->nlines());
        quad_range_weight_(i).resize(model.parameters->nquads());

        intensities_(i).resize(length, width);
        delta_tau_(i).resize(length, width);
        S_curr_(i).resize(length, width);
        S_next_(i).resize(length, width);

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

    //Finally also trace the rays in advance to prune the unnecessary ones.
    get_static_rays_to_trace(model);

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
    const Size      o,//origin point of ray
    const Size      rdir,//ray direction to trace
    const Size      rsav,//ray index save direction
    const Size      rayidx)//for indexing the rays (counts from 0 to total number of rays traced through the domain
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
        Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);

        if (dist2<min_ray_distsqr(rsav,nxt))
        {
            min_ray_distsqr(rsav,nxt)=dist2;
            closest_ray(rsav,nxt)=rayidx;
        }

        Size crt = o;

        while (geometry.not_on_boundary(nxt))
        {
            crt =       nxt;
            nxt = geometry.get_next(o, rdir, nxt, Z, dZ);

            n_rays_through_point(rsav,nxt)++;

            // get distance and check if closest ray
            Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);

            if (dist2<min_ray_distsqr(rsav,nxt))
            {
                min_ray_distsqr(rsav,nxt)=dist2;
                closest_ray(rsav,nxt)=rayidx;
            }
        }
    }
}

// For all direction, determines a ray covering of the points
inline void Solver :: get_static_rays_to_trace (Model& model)
{
    accelerated_for (rr, model.parameters->hnrays(),
    {
        Size n_rays_to_trace=0;
        const Size ar=model.geometry.rays.antipod[rr];
        for (Size o=0; o<model.parameters->npoints(); o++)
        {
            // std::cout<<"n_rays_through_point: "<<n_rays_through_point(rr,o)<<std::endl;
            // TODO?: maybe replace with some boolean-like thing, as we do not actually care (aside from debugging)
            // how much rays are traced through a point
            if (n_rays_through_point(rr,o)>0)
            {continue;}

            //seemingly no ray has been traced through this point, so we must trace a ray through it
            points_to_trace_ray_through[rr][n_rays_to_trace]=o;
            //tracing rays is symmetric, so only keep for the first half of the ray directions

            //trace ray through point
            n_rays_through_point(rr,o)++;
            //antipod has exactly same number of rays through the point, so do not save

            // For generality, I assign a ray index to each ray of a given direction rsav
            // Also the ray direction is identified with the forward dir
            closest_ray(rr,o)=rr;
            min_ray_distsqr(rr,o)=0.0;

            //now trace rest of rays
            trace_ray_points(model.geometry, o, rr, rr, n_rays_to_trace);
            trace_ray_points(model.geometry, o, ar, rr, n_rays_to_trace);

            n_rays_to_trace++;
        }
        // n_points_to_trace_ray_through[rr]=n_rays_to_trace;
        points_to_trace_ray_through[rr].resize(n_rays_to_trace);//and now the correct size, instead of parameters->npoints()

    })

    //debug print stuff
    for (Size rr=0; rr<model.parameters->hnrays(); rr++)
    {
        std::cout<<"rr: "<<rr<<" size points_to_trace_ray_through: "<<points_to_trace_ray_through[rr].size()<<std::endl;
        for (Size idx=0; idx<points_to_trace_ray_through[rr].size(); idx++)
        {
            // std::cout<<"point: "<<points_to_trace_ray_through[rr][idx]<<std::endl;
        }
        // std::cout<<"number of rays per point"<<std::endl;
        for (Size p=0; p<model.parameters->npoints(); p++)
        {
            // std::cout<<"point: "<<p<<"#: "<<n_rays_through_point(rr, p)<<std::endl;
        }

    }

}


// // Splits the frequencies into different parts, such that the different parts are sufficiently seperated
// // TODO: In single line approx, this step should be ignored, as all lines are far enough from eachother
// //ASSUMES FREQS SORTED
// // p should be next point on the ray
// //NOTE TO SELF: I ASSUMED FOR SIMPLICITY THAT I SPLIT IT A SINGLE RAY AT A TIME; NOT DIFFERENTLY FOR EACH POINT ON THE RAY
// //TODO FIXME: USE THE RAY TO SPLIT THE FREQUENCIES IN A GLOBAL BOUND MANNER
// // as a consequence, handling frequency quadratures very close to eachother should be automatically included
// //Note: this can mean however that the frequency derivative terms might not be computed that accurately if the lines separate
// inline void Solver :: split_frequencies (Model& model, Size p)
// {
//     freqsplits_()[0]=0;
//     n_freqsplits_()=1;
//     Real currfreq=model.radiation.frequencies.nu(p, 0);
//     Real deltafreq=0.0;
//     //TODO GET MIN LINE WIDTH
//     for (freqid=1; freqid<model.parameters->nfreqs(), freqid++)
//     {
//         deltafreq=model.radiation.frequencies.nu(p, freqid)-currfreq;
//         if (deltafreq>FACTORTODO*LINEWIDTHTODO)
//         {
//             freqsplits_()[n_freqsplits_()]=freqid;
//             n_freqsplits_()++;
//         }
//         currfreq=model.radiation.frequencies.nu(p, freqid);
//     }
//     freqsplits_()[n_freqsplits_()]=model.parameters->nfreqs();//for simplicity, also add the last split.. makes it easier to get rightmost value in a split
// }


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
        //Starting from the lowest frequency
        Size curr_freq_idx = 0;
        Size next_freq_idx = 0;
        const Size max_freq_idx = model.parameters->nfreqs()-1;
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        while (curr_freq_idx<=max_freq_idx && next_freq_idx<=max_freq_idx)
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
    Vector<unsigned char>& line_quad_discdir_overlap=line_quad_discdir_overlap_();
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
        while (quadidx<model.parameters->nquads()-2)
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }

        while (quadidx<model.parameters->nquads())
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

        while (quadidx<model.parameters->nquads())
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }
    }

    //Dummy initialization, as we do not know a priori what the first counted quadrature will be
    Real leftbound=0.0;
    Real rightbound=0.0;
    Size leftboundidx=0;
    Size rightboundidx=0;

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
inline void Solver :: set_implicit_boundary_frequencies(Model& model, const Size nextpoint, const Size nextpointonrayindex, Real shift_next,
                std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index, bool is_upward_disc)
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
    Vector<Real>& left_bound=left_bound_();
    Vector<Real>& right_bound=right_bound_();
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
            if (quadidx>=model.parameters->nquads()-2&&!line_quad_overlap[lineidx])
            {
                //then is required boundary condition
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
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
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
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
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
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
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
            }
        }
    }
}

// Matches the overlapping boundary conditions, using the currently computed ranges
inline void Solver :: match_overlapping_boundary_conditions(Model& model, const Size currpoint, const Size curr_point_on_ray_index,
                                                            const Real curr_shift, std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index)
{
    //just iterate simultaneously over the ranges and frequencies
    //TPDO first define the indices corresponding to the ranges//DONE
    Vector<Real>& left_bound=left_bound_();//Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound=right_bound_();//Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges=nb_ranges_();//contains the number of ranges
    Vector<Size>& left_bound_index=left_bound_index_();//Specifies the corresponding freq index to the left bounds of the ranges
    Vector<Size>& right_bound_index=right_bound_index_();//Specifies the corresponding freq index to the right bounds of the ranges

    Size curr_range_index=0;
    Size curr_freq_index=left_bound[curr_range_index];//contains the current frequency index
    //Assumption: at least a single line exists -> a single range exists
    Real curr_range_min=left_bound[curr_range_index];
    Real curr_range_max=right_bound[curr_range_index];
    Real curr_range_freq=curr_range_min;


    //TODO: also add frequencies!
    //sorted iterating over all boundary frequencies
    //TODO: just get iterator over multimap and simulataneously check for end of ranges too!
    std::multimap<Real, std::tuple<Size, Size>>::iterator it=multimap_freq_to_bdy_index.begin();
    while (it!=multimap_freq_to_bdy_index.end()&&curr_range_index<nb_ranges)
    {
        const Real bdy_freq=it->first;//get boundary freq
        const Size bdy_point_on_ray_idx=std::get<0>(it->second);//get corresponding point index on ray
        const Size bdy_freq_idx=std::get<1>(it->second);//get corresponding freq index

        //compare the frequency to the max freq of current range; if larger, then we need to compare against the next range
        while (right_bound[curr_range_index]<bdy_freq)
        {
            if (curr_range_index==nb_ranges)
            {
                //if curr_range_index==nb_ranges, then all remaining boundary frequencies lie outside of the ranges
                //so we do not need to waste any more time on trying to match them
                return;
            }
            curr_range_index++;
            //update info about ranges
            curr_range_min=left_bound[curr_range_index];
            curr_range_max=right_bound[curr_range_index];
            curr_range_freq=curr_range_min;
            curr_freq_index=left_bound_index[curr_range_index];//stores the freq index of curr point
        }
        //now it is guaranteed that rightbound>=bdy_freq, so we only need to check whether leftbound<=bdy_freq
        if (curr_range_min<=bdy_freq)
        {
            //Do the boundary stuff
            //but first check which freqs are exactly in front or behind the bdy freq
            //Start with the leftmost index of the range

            //while(&& <) looping
            //Find which freq bounds at curr_point exactly correspond to the bdy freq (enclosing it)
            while (bdy_freq>curr_range_freq&&curr_range_freq<curr_range_max)
            {
                curr_freq_index++;
                Real curr_range_freq=model.radiation.frequencies.nu(currpoint, curr_freq_index)*curr_shift;
            }
            Size left_curr_freq_idx=curr_freq_index-1;
            Size right_curr_freq_idx=curr_freq_index;
            //the curr range freq should now be larger/equal to the bdy freq
            //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as index for interpolation, so use curr_index as left bound instead
            if (curr_freq_index==left_bound_index[curr_range_index])
            {
                left_curr_freq_idx=curr_freq_index;
                right_curr_freq_idx=curr_freq_index+1;
            }
            const Real deltafreq=(model.radiation.frequencies.nu(currpoint, right_curr_freq_idx)
                                  -model.radiation.frequencies.nu(currpoint, left_curr_freq_idx))*curr_shift;

            //Set starting index correctly
            start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[0]=curr_point_on_ray_index;
            start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1]=left_curr_freq_idx;

            //Set frequency derivative correctly

            //Set explicit coefficents to TWICE the normal 1/Δν value (as we are only treating the explicit part)
            dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=-2/deltafreq;
            dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=+2/deltafreq;
            dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

            dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
            dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;
            dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;
            //as the coefficient for the last on is 0, no effort should be made to renumber that last index

            //Implicit part coefs are set to 0. (as long as the corresponding indicies are in bounds, i do not particularly care about the exact values)
            dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
            dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
            dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
            //exact indices do not matter, as long as they are in bounds
            dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
            dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
            dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;

            //set dtau to approx 0
            delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx)=COMOVING_MIN_DTAU;//should be small enough, but not small enough to crash my solver (due to /Δτ^2 necessary)
            std::cout<<"bdy delta_tau: "<<COMOVING_MIN_DTAU<<std::endl;
            //set S? (nah just ignore the existence, as dtau≃0)

            //pop value from iterator
            it=multimap_freq_to_bdy_index.erase(it);
            //TODO NOT HERE: for the INITIAL BOUNDARY CONDITIONS, maybe let the function choose the TYPE OF BOUNDARY CONDITION to use!!!
        }
        else
        {
            it++;
        }
    }
}

//As currpoint lies on the boundary, use currpoint to derive the necessary boundary conditions!
inline void Solver :: set_initial_boundary_conditions(Model& model, const Size currpoint, const Real curr_shift, std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index)
{
    std::multimap<Real, std::tuple<Size, Size>>::iterator it=multimap_freq_to_bdy_index.begin();
    while (it!=multimap_freq_to_bdy_index.end())
    {
        const Real bdy_freq=it->first;//get boundary freq
        const Size bdy_point_on_ray_idx=std::get<0>(it->second);//get corresponding point index on ray
        const Size bdy_freq_idx=std::get<1>(it->second);//get corresponding freq index

        //For using somewhat reasonable frequency indices
        Size curr_freq_idx=start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1];

        //Compute boundary intensity
        Real bdy_intensity=boundary_intensity(model, currpoint, bdy_freq);

        //Set the frequency derivative coefficients to 0
        dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

        dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;
        dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;
        dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;

        dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

        dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;

        // Set S and dtau such that the resulting intensity will be (approximately) equal to the initial boundary intensity
        delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx)=50.0;//should be large enough
        S_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_intensity;
        S_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_intensity;
        std::cout<<"bdy intensity: "<<bdy_intensity<<std::endl;

        //pop value from iterator
        it=multimap_freq_to_bdy_index.erase(it);
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


//Computes the setup and the boundary conditions for the comoving solver for the forward ray direction
//Assumes we have already traced the ray
//TODO: do MORE specialization for the single line approx (taking inspiration from far too complex deque stuff?)
template<ApproximationType approx>
inline void Solver :: comoving_ray_bdy_setup_forward(Model& model)
{
    std::cout<<"entering forward setup"<<std::endl;
    //Note: this assumes forward ray, not backward ray TODO IMPLEMENT OTHER DIRECTION
    Matrix<Real>& intensities=intensities_();

    std::cout<<"after renaming intensities_"<<std::endl;
    //TODO set o=index of first point on ray
    Size first_index=first_();//index of first point on ray
    Size shift_first=2.0-shift_()[first_index];
    Size rayposidx=last_();//ray position index -> point index through nr_()[rayposidx]
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        intensities(first_index, freqid)=boundary_intensity(model, nr_()[first_index], model.radiation.frequencies.nu(nr_()[first_index], freqid)*shift_first);
    }

    std::cout<<"after setting first intensities"<<std::endl;

    //from first_+1 to last_ index, trace the ray in reverse direction to treat the boundary conditions

    //dummy initialization; just need some references to Reals for get_eta_and_chi
    //Previously defined vectors are not of the right size for me (O(npointsonray) instead of O(nfreqs))//
    //Wrong: eta_c_, ... are the correct dimension, however I see no reason to use them. (as the values will be directly used)
    //So I just compute per point individually (however the compiler should notice and loop unrolling could still happen)
    Real eta_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real chi_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real eta_curr=0.0;
    Real chi_curr=0.0;

    //helper values for computing the freq derivative terms
    Real dfreqsmall=0.0;
    Real dfreqlarge=0.0;

    //POSSIBLE OPTIMIZATION: only query them once, moving ..._curr to ..._next after each rayposidx change
    //But then this should change to a vector holding these values...

    std::multimap<Real, std::tuple<Size, Size>> multimap_freq_to_bdy_index;//boundary conditions can overlap, so multimap it is


    std::cout<<"starting iteration"<<std::endl;
    //Yes, for determining boundary conditions, we are going backwards; this seems strange, but is a result of easier treatment of the spatial boundary conditions at the end
    while (rayposidx>first_())
    {
        const Size nextpoint=nr_()[rayposidx];
        const Size currpoint=nr_()[rayposidx-1];
        const Real dZ=dZ_()[rayposidx-1];
        //TODO: for less mistakes, redefine rayposidx, rayposidx-1 using different names

        const Real shift_next=2.0-shift_()[rayposidx];//shift is by default defined backwards for the forward ray direction; so under the non-relativistic approx, we can just do 2-shift to get the shift in the other direction
        const Real shift_curr=2.0-shift_()[rayposidx-1];//TODO: IS DIFFERENT FOR BACKWARD RAY

        //For the sake of accurately putting boundary conditions, this should point in the correct direction; otherwise we might accidentally put boundary conditions somewhat closer to the middle of a line
        const bool is_upward_disc=(shift_next>=shift_curr);

        //TODO: maybe refactor loop from this point onwards; should be the same for both forward and backward rays

        std::cout<<"before matching freq indices"<<std::endl;

        //somewhere, we should also compute all of S, Δτ, dI/dν stuff
        //this first part might be appropriate, as we can as well overwrite all this stuff later for the boundary conditions
        //err, Δτ is not available until we know which freq matches which other freq
        //First, we match all frequency indices as good as possible.
        match_frequency_indices(model, nextpoint, currpoint, rayposidx, rayposidx-1, is_upward_disc);//O(nfreqs)

        std::cout<<"before setting default computation values"<<std::endl;

        //Now compute all default stuff, computing bogus for the non-matched indices (but as these correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Real nextfreq=model.radiation.frequencies.nu(nextpoint, next_freq_idx);//=comoving frame freq
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq=model.radiation.frequencies.nu(currpoint, curr_freq_idx);//
            //technically, we also need to read the point index from this start_indices_ (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to currpoint at this moment in time
            const Size nextlineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, next_freq_idx);
            get_eta_and_chi<approx>(model, nextpoint, nextlineidx, nextfreq, eta_next, chi_next);//no shift necessary, as the frequencies are matched
            //Possible optimization, only query eta/chi_curr once for each point, as the next value is the same either way when traversing through this ray
            const Size currlineidx=model.radiation.frequencies.corresponding_line_matrix(currpoint, curr_freq_idx);
            get_eta_and_chi<approx>(model, currpoint, currlineidx, currfreq, eta_curr, chi_curr);
            const Real Snext=eta_next/chi_next;
            const Real Scurr=eta_curr/chi_curr;
            S_next_()(rayposidx, next_freq_idx)=Snext;
            S_curr_()(rayposidx, next_freq_idx)=Scurr;
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            const Real dtau = std::max(trap(chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
            delta_tau_()(rayposidx, next_freq_idx)=dtau;
            std::cout<<"setting dtau: "<<delta_tau_()(rayposidx, next_freq_idx)<<std::endl;
            std::cout<<"chi_next: "<<chi_next<<std::endl;
            std::cout<<"chi_curr: "<<chi_curr<<std::endl;
            std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        }

        std::cout<<"setting derivative coeffs"<<std::endl;

        //This part precomputes the frequency derivative coefficients (and corresponding indices), conveniently ignoring the fact that for some outermost line quadrature frequencies, we must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        //The discretization is evidently depends on freq disc direction
        if (is_upward_disc)
        {
            std::cout<<"is_upward_disc"<<std::endl;
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            //FIXME: replace with while loop to deal with nfreqs < 3!
            for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs()-2; next_freq_idx++)
            {
                std::cout<<"get start indices"<<std::endl;
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                std::cout<<"getting start_indices_"<<std::endl;
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices
                Size curr_freq_idxp1=(curr_freq_idx+1)%model.parameters->nfreqs();//index+1
                Size curr_freq_idxp2=(curr_freq_idx+2)%model.parameters->nfreqs();//index+2
                std::cout<<"freq+2: "<<curr_freq_idxp2<<" nfreqs: "<<model.parameters->nfreqs()<<std::endl;
                std::cout<<"computing curr coeffs"<<std::endl;
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                dfreqsmall=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
                            *shift_curr;
                dfreqlarge=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
                            *shift_curr;
                std::cout<<"setting curr coeffs and indices"<<std::endl;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                std::cout<<"computing next coeffs"<<std::endl;

                //And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx+1))
                            *shift_next;
                dfreqlarge=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx))
                            *shift_next;

                std::cout<<"setting next coeffs and indices"<<std::endl;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                std::cout<<"setting next indices"<<std::endl;
                std::cout<<"rayposidx: "<<rayposidx<<std::endl;
                std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx+2;
                std::cout<<"flag 3"<<std::endl;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx+1;
                std::cout<<"flag 2"<<std::endl;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;
                std::cout<<"flag 1"<<std::endl;

            }
            //and also just set some irrelevant values for the final boundary points, as these should get overwritten?
            //Wait do we then even need to set it here?
            //TODO: figure out whether I need to set some random values here!
        }
        else
        {
            std::cout<<"!is_upward_disc"<<std::endl;
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            //Assumes nfreqs>=3
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
                *shift_curr;
                dfreqlarge=-(model.radiation.frequencies.nu(currpoint, curr_freq_idx)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;
                dfreqlarge=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }

        //Note to self: we can just use all default freq derivative stuff, except for the 2 outermost points on the grid.
        // Because I define every computation in terms of the next point index, I will have no problems with non-existant/too far frequencies, as I have made the computed curr regions slightly smaller
        //Note: it might seem we will try to use oob values, but these should be overwritten anyway
        //CHECK THIS FOR A SIMPLE EXAMPLE!!!

        std::cout<<"getting line ranges"<<std::endl;

        //now start classifying the new boundary points
        //For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc, shift_curr);//for the current point, we definitely need the exact ranges (fiddled with to ensure enough bdy conditions)

        std::cout<<"getting overlapping lines"<<std::endl;
        // and compute the overlap between lines on next_point
        get_overlapping_lines(model, nextpoint, is_upward_disc);//technically, we should compute which lines overlap only in this part

        std::cout<<"setting implicit boundary freqs"<<std::endl;
        //First, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be treated as boundary points
        set_implicit_boundary_frequencies(model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);

        std::cout<<"matching implicit boundary freqs"<<std::endl;
        // And now we check what boundary conditions need to be evaluated using the intensities at curr_point
        match_overlapping_boundary_conditions(model, currpoint, rayposidx-1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx--;
    }
    std::cout<<"setting initial bounds"<<std::endl;
    //After going through all points, the remaining boundary frequencies have not matched any frequencies near any line,
    // so we will use initial boundary conditions computed at the first point to put on the ray
    //note: currpoint should be point on boundary of domain; in this way, we can simply add the corresponding boundary conditions
    set_initial_boundary_conditions(model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
}


//Computes the setup and the boundary conditions for the comoving solver for the backward ray direction
//Assumes we have already traced the ray
//TODO: do MORE specialization for the single line approx (taking inspiration from far too complex deque stuff?)
template<ApproximationType approx>
inline void Solver :: comoving_ray_bdy_setup_backward(Model& model)
{
    //Note: this assumes forward ray, not backward ray TODO IMPLEMENT OTHER DIRECTION
    Matrix<Real>& intensities=intensities_();
    //TODO set o=index of first point on ray
    Size first_index=last_();//index of first point on ray
    Size shift_first=shift_()[first_index];
    Size rayposidx=first_();//ray position index -> point index through nr_()[rayposidx]
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        intensities(first_index, freqid)=boundary_intensity(model, nr_()[first_index], model.radiation.frequencies.nu(nr_()[first_index], freqid)*shift_first);
    }

    //from first_+1 to last_ index, trace the ray in reverse direction to treat the boundary conditions

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
    //But then this should change to a vector holding these values...

    std::multimap<Real, std::tuple<Size, Size>> multimap_freq_to_bdy_index;//boundary conditions can overlap, so multimap it is
    //Yes, for determining boundary conditions, we are going backwards; this seems strange, but is a result of easier treatment of the spatial boundary conditions at the end
    while (rayposidx<last_())
    {
        const Size nextpoint=nr_()[rayposidx];
        const Size currpoint=nr_()[rayposidx+1];
        const Real dZ=dZ_()[rayposidx];//dZ is stored somewhat finnicky, in locations [first_(),last_()[ (so not including last_())
        //TODO: for less mistakes, redefine rayposidx, rayposidx+1 using different names

        const Real shift_next=shift_()[rayposidx];//shift is by default defined backwards for the forward ray direction; so under the non-relativistic approx, we can just do 2-shift to get the shift in the other direction
        const Real shift_curr=shift_()[rayposidx+1];//TODO: IS DIFFERENT FOR BACKWARD RAY

        //For the sake of accurately putting boundary conditions, this should point in the correct direction; otherwise we might accidentally put boundary conditions somewhat closer to the middle of a line
        const bool is_upward_disc=(shift_next>=shift_curr);

        //TODO: maybe refactor loop from this point onwards; should be the same for both forward and backward rays

        //somewhere, we should also compute all of S, Δτ, dI/dν stuff
        //this first part might be appropriate, as we can as well overwrite all this stuff later for the boundary conditions
        //err, Δτ is not available until we know which freq matches which other freq
        //First, we match all frequency indices as good as possible.
        match_frequency_indices(model, nextpoint, currpoint, rayposidx, rayposidx+1, is_upward_disc);//O(nfreqs)

        //Now compute all default stuff, computing bogus for the non-matched indices (but as these correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Real nextfreq=model.radiation.frequencies.nu(nextpoint, next_freq_idx);//=comoving frame freq
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq=model.radiation.frequencies.nu(currpoint, curr_freq_idx);//
            //technically, we also need to read the point index from this start_indices_ (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to currpoint at this moment in time
            const Size nextlineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, next_freq_idx);
            get_eta_and_chi<approx>(model, nextpoint, nextlineidx, nextfreq, eta_next, chi_next);//no shift necessary, as the frequencies are matched
            //Possible optimization, only query eta/chi_curr once for each point, as the next value is the same either way when traversing through this ray
            const Size currlineidx=model.radiation.frequencies.corresponding_line_matrix(currpoint, curr_freq_idx);
            get_eta_and_chi<approx>(model, currpoint, currlineidx, currfreq, eta_curr, chi_curr);
            const Real Snext=eta_next/chi_next;
            const Real Scurr=eta_curr/chi_curr;
            S_next_()(rayposidx, next_freq_idx)=Snext;
            S_curr_()(rayposidx, next_freq_idx)=Scurr;
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            const Real dtau = std::max(trap (chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
            delta_tau_()(rayposidx, next_freq_idx)=dtau;
            std::cout<<"setting dtau: "<<delta_tau_()(rayposidx, next_freq_idx)<<std::endl;
            std::cout<<"chi_next: "<<chi_next<<std::endl;
            std::cout<<"chi_curr: "<<chi_curr<<std::endl;
            std::cout<<"dZ: "<<dZ<<std::endl;
            std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        }

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
                            *shift_curr;
                dfreqlarge=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
                            *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx+1))
                            *shift_next;
                dfreqlarge=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx))
                            *shift_next;

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
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side (switched order of term, making sign clearer; just compare to upward disc version)
                dfreqsmall=-(model.radiation.frequencies.nu(currpoint, curr_freq_idxm1)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dfreqlarge=-(model.radiation.frequencies.nu(currpoint, curr_freq_idx)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;
                dfreqlarge=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }

        //Note to self: we can just use all default freq derivative stuff, except for the 2 outermost points on the grid.
        // Because I define every computation in terms of the next point index, I will have no problems with non-existant/too far frequencies, as I have made the computed curr regions slightly smaller
        //Note: it might seem we will try to use oob values, but these should be overwritten anyway
        //CHECK THIS FOR A SIMPLE EXAMPLE!!!

        //now start classifying the new boundary points
        //For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc, shift_curr);//for the current point, we definitely need the exact ranges (fiddled with to ensure enough bdy conditions)
        // and compute the overlap between lines on next_point
        get_overlapping_lines(model, nextpoint, is_upward_disc);//technically, we should compute which lines overlap only in this part

        //First, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be treated as boundary points
        set_implicit_boundary_frequencies(model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);
        // And now we check what boundary conditions need to be evaluated using the intensities at curr_point
        match_overlapping_boundary_conditions(model, currpoint, rayposidx+1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx++;
    }
    //After going through all points, the remaining boundary frequencies have not matched any frequencies near any line,
    // so we will use initial boundary conditions computed at the first point to put on the ray
    //note: nr_()[first_index] should be point on boundary of domain; in this way, we can simply add the corresponding boundary conditions
    set_initial_boundary_conditions(model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
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

    std::cout<<"hnrays: "<<model.parameters->hnrays()<<std::endl;
    std::cout<<"number threads: "<<paracabs::multi_threading::n_threads ()<<std::endl;
    // For each ray, solve the radiative transfer equation in a comoving manner
    accelerated_for (rr, model.parameters->hnrays(),
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        // const Real     wt = model.geometry.rays.weight   [rr] * two;
        // const Vector3D nn = model.geometry.rays.direction[rr];

        std::cout << "--- rr = " << rr << std::endl;
        std::cout<<"number threads: "<<paracabs::multi_threading::n_threads ()<<std::endl;

        std::cout<<"thread id: "<<paracabs::multi_threading::thread_id()<<std::endl;

        //for every ray to trace
        const Size n_rays_to_trace=points_to_trace_ray_through.size();
        for (Size rayidx=0; rayidx<n_rays_to_trace; rayidx++)
        {
        // accelerated_for (rayidx, n_rays_to_trace,
        // {
            //get corresponding origins of the rays
            const Size o = points_to_trace_ray_through[rr][rayidx];
            const double dshift_max = get_dshift_max (model, o);
            std::cout<<"here"<<std::endl;
            //FIXME: using dshift_max computed from a single point only works when the line width does not change too much; however, all other solvers are already using this bad approximation
            //trace and solve over the ray in both directions, incrementing J as necessary
            solve_comoving_order_2_sparse<approx>(model, o, rr, rayidx, dshift_max);//solves both for forward and backward ray
            //complicated to decouple ray tracing from solving, so more logic is moved inside this functions
        // })
        }
    })

    // // TODO REPLACE/DELETE THIS?
    // // For each ray, solve transfer equation
    // for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    // {
    //     const Size     ar = model.geometry.rays.antipod  [rr];
    //     const Real     wt = model.geometry.rays.weight   [rr] * two;
    //     const Vector3D nn = model.geometry.rays.direction[rr];
    //
    //     cout << "--- rr = " << rr << endl;
    //
    //     for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    //     {
    //         threaded_for (o, model.parameters->npoints(),
    //         {
    //             const Real dshift_max = get_dshift_max (model, o);
    //
    //             nr_   ()[centre] = o;
    //             shift_()[centre] = 1.0;
    //
    //             first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    //             last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
    //             n_tot_() = (last_()+1) - first_();
    //
    //             if (n_tot_() > 1)
    //             {
    //                 for (Size k = 0; k < lspec.linedata.nrad; k++)
    //                 {
    //                     // Integrate over the line
    //                     for (Size z = 0; z < model.parameters->nquads(); z++)
    //                     {
    //                         solve_feautrier_order_2 <approx> (model, o,  lspec.nr_line[o][k][z]);
    //
    //                         lspec.J(o,k) += lspec.quadrature.weights[z] * wt * Su_()[centre];
    //
    //                         update_Lambda           (model, rr, lspec.nr_line[o][k][z]);
    //                     }
    //                 }
    //             }
    //             else
    //             {
    //                 for (Size k = 0; k < lspec.linedata.nrad; k++)
    //                 {
    //                     // Integrate over the line
    //                     for (Size z = 0; z < model.parameters->nquads(); z++)
    //                     {
    //                         lspec.J(o,k) += lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
    //                     }
    //                 }
    //             }
    //         })
    //     }
    // }
}


//As stepping a single step should be exactly the same in both directions on the ray, we might as well refactor it out
// //rayposidx==rayposidx_curr_point+-1
//Ergo rayposidx stand for the ray position index of the next point
//rr direction index necessary for determining whether to add intensity to J//TODO? replace with bool denoting whether to add it (precompute if clause somewhere else)
//TODO: figure out whether rayposidx_currpoint is needed!! As it should be replaced with start_indices_...
inline void Solver :: solve_comoving_single_step (Model& model, const Size rayposidx, const Size rayidx, const Size rr, const bool is_upward_disc, const bool forward_ray)
{
    Vector<Size>& nr=nr_();//stores the exact point indices
    Vector<unsigned char>& real_pt=real_pt_();//stores whether each point is real point (not added extra for interpolation purposes)
    Vector<double>& shift=shift_();//stores the shifts versus the static frame; warning: for forward rays, this should be changed to 2.0-shift; for backward rays, it is correct
    Matrix<Vector<Size>>& start_indices=start_indices_();//for getting the previous indices associated with computing the intensity for the next point
    Matrix<Real>& delta_tau=delta_tau_();//stores the optical depths
    Matrix<Real>& S_curr=S_curr_();
    Matrix<Real>& S_next=S_next_();
    Matrix<Real>& intensities=intensities_();

    Matrix<Real>& dIdnu_coef1_curr=dIdnu_coef1_curr_();//coefficient of the freq point itself
    Matrix<Real>& dIdnu_coef2_curr=dIdnu_coef2_curr_();
    Matrix<Real>& dIdnu_coef3_curr=dIdnu_coef3_curr_();

    Matrix<Size>& dIdnu_index1_curr=dIdnu_index1_curr_();
    Matrix<Size>& dIdnu_index2_curr=dIdnu_index2_curr_();
    Matrix<Size>& dIdnu_index3_curr=dIdnu_index3_curr_();

    Matrix<Real>& dIdnu_coef1_next=dIdnu_coef1_next_();
    Matrix<Real>& dIdnu_coef2_next=dIdnu_coef2_next_();
    Matrix<Real>& dIdnu_coef3_next=dIdnu_coef3_next_();

    Matrix<Size>& dIdnu_index1_next=dIdnu_index1_next_();
    Matrix<Size>& dIdnu_index2_next=dIdnu_index2_next_();
    Matrix<Size>& dIdnu_index3_next=dIdnu_index3_next_();

    std::cout<<"after renaming all vars"<<std::endl;
    std::cout<<"rayposidx=next point on ray idx= "<<rayposidx<<std::endl;

    //rayposidx represents the next position on the ray
    const Size nextpointidx=nr[rayposidx];
    const Real next_shift=(forward_ray) ? 2.0-shift[rayposidx] : shift[rayposidx];//absolute shifts, versus static frame!

    std::cout<<"after getting shift"<<std::endl;

    //Do the explicit part
    for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
    {
        std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
        //Get the point indices
        const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
        std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
        const Size curr_point_idx=nr[curr_point_on_ray_index];
        const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
        const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame

        std::cout<<"after computing shift"<<std::endl;

        const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
        const Real dtau=delta_tau(rayposidx, next_freq_idx);
        std::cout<<"dtau: "<<dtau<<std::endl;
        std::cout<<"rayposidx: "<<rayposidx<<std::endl;

        std::cout<<"after computing dtau"<<std::endl;

        const Real expl_term=(-expm1(-dtau)-dtau*exp(-dtau))/dtau;
        std::cout<<"expl term: "<<expl_term<<std::endl;
        std::cout<<"curr source: "<<S_curr(rayposidx, next_freq_idx)<<std::endl;

        const Real expl_freq_der=dIdnu_coef1_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))
        +dIdnu_coef2_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))
        +dIdnu_coef3_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx));
        std::cout<<"expl_freq_der: "<<expl_freq_der<<std::endl;
        //TODO: actually do the explicit part
        intensities(rayposidx, next_freq_idx)=intensities(curr_point_on_ray_index, curr_freq_idx)*exp(-dtau)
        +expl_term*S_curr(rayposidx, next_freq_idx) //source term
        +expl_term*expl_freq_der*deltanu/dtau;//frequency derivative term
        std::cout<<"explicit part intensities: "<<intensities(rayposidx, next_freq_idx)<<std::endl;
    }

    std::cout<<"after explicit part"<<std::endl;

    //Implicit part ordering depends on the discretization direction, use whether we have an upward discretization instead
    // const Real rel_doppler_shift=shift[rayposidx]-shift[rayposidx+1];

    // if (rel_doppler_shift>=0)
    if (is_upward_disc)
    {
        std::cout<<"is_upward_disc"<<std::endl;
        //in case of positive doppler shifts, we need to start from the largest frequencies (where we have extra boundary conditions on the end)
        //reverse loop over unsigned int, so indices +1
        for (Size next_freq_idx=model.parameters->nfreqs(); next_freq_idx>0; next_freq_idx--)
        {
            //TODO: actually do the implicit part
            const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx-1)[0];
            const Size curr_point_idx=nr[curr_point_on_ray_index];
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx-1)[1];
            const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx-1)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real dtau=delta_tau(rayposidx, next_freq_idx-1);
            std::cout<<"dtau: "<<dtau<<std::endl;
            std::cout<<"next freq idx: "<<next_freq_idx-1<<std::endl;
            std::cout<<"rayposidx: "<<rayposidx<<std::endl;
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            std::cout<<"impl term: "<<impl_term<<std::endl;
            std::cout<<"next source: "<<S_next(rayposidx, next_freq_idx-1)<<std::endl;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=dIdnu_coef2_next(rayposidx, next_freq_idx-1)*intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx-1))
            +dIdnu_coef3_next(rayposidx, next_freq_idx-1)*intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx-1));
            std::cout<<"impl freq der: "<<impl_freq_der<<std::endl;

            intensities(rayposidx, next_freq_idx-1)=(intensities(rayposidx, next_freq_idx-1)
            +impl_term*S_next(rayposidx, next_freq_idx-1)//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-dIdnu_coef1_next(rayposidx, next_freq_idx-1)*deltanu/dtau); //freq derivative term
        }
    }
    else
    {
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            //TODO: actually do the implicit part
            const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
            const Size curr_point_idx=nr[curr_point_on_ray_index];
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real dtau=delta_tau(rayposidx, next_freq_idx);
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=dIdnu_coef2_next(rayposidx, next_freq_idx)*intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx))
            +dIdnu_coef3_next(rayposidx, next_freq_idx)*intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx));

            intensities(rayposidx, next_freq_idx)=(intensities(rayposidx, next_freq_idx)
            +impl_term*S_next(rayposidx, next_freq_idx)//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-dIdnu_coef1_next(rayposidx, next_freq_idx)*deltanu/dtau); //freq derivative term
        }
    }

    std::cout<<"incrementing J"<<std::endl;
    std::cout<<"real_pt[rayposidx]: "<<real_pt[rayposidx]<<std::endl;
    if (real_pt[rayposidx])
    {
        std::cout<<"is real point"<<std::endl;
    }else{
        std::cout<<"is not real point"<<std::endl;
    }
    std::cout<<"closest ray(...): "<<closest_ray(rr, nextpointidx)<<std::endl;
    //? if condition always says no!

    //Finally increment J if and real point check if closest ray
    if (real_pt[rayposidx]&&closest_ray(rr, nextpointidx)==rayidx)
    {
        std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        std::cout<<"adding to J at rayposidx: "<<rayposidx<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< quadrature number corresponding to frequency
            // TODO: replace with matrix variant when replacing these Vectors
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            std::cout<<"freqid: "<<freqid<<std::endl;
            std::cout<<"adding: "<<lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid)<<std::endl;
            std::cout<<"intensity: "<<intensities(rayposidx, freqid)<<std::endl;
            lspec.J(nextpointidx,k) += lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid);// Su_()[centre];
            //TODO: compute lambda term
        }
    }
}


template<ApproximationType approx>
inline void Solver :: solve_comoving_order_2_sparse (
      Model& model,
      const Size o,//ray origin point
      const Size r,//ray direction index
      const Size rayidx,//ray trace index
      const double dshift_max)
{
    const Size rr=r;
    const Size ar=model.geometry.rays.antipod[rr];
    //FIXME: during setup, check if not spherically symmetric model!!

    // Vector<Real>& eta_c = eta_c_();//current emissivity for all freqs
    // Vector<Real>& eta_n = eta_n_();//next emissivity    for all freqs
    //
    // Vector<Real>& chi_c = chi_c_();//current opacity    for all freqs
    // Vector<Real>& chi_n = chi_n_();//next opacity       for all freqs

    // Vector<Real>& curr_intensity = curr_intensity_();//current intensity for all freqs
    // Vector<Real>& next_intensity = next_intensity_();//storage for next intensity for all freqs

    Matrix<Real>& intensities=intensities_();//will contain all computed intensities

    Vector<unsigned char>& real_pt= real_pt_();//for denoting whether the given point is a real point
    Vector<Size>& nr=nr_();

    Vector<double>& shift = shift_();// contains the doppler shifts for all points on the ray
    //for the forward ray, do 2.0-shift to get the shift direction I want for my frequencies
    //for the backward ray, this should be correct
    // TODO: fill in yourself

    // const Real     wt = model.geometry.rays.weight   [rr];//ray direction weight

    double Z = 0.0; //distance along ray
    double dZ= 0.0; //last distance increment

    std::cout<<"before tracing rays"<<std::endl;

    //trace the ray, getting all points on the ray and indicating whether they are real
    first_() = trace_ray_indicate_point(model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    last_ () = trace_ray_indicate_point(model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
    //also set ray start as real point; accidentally forgot
    real_pt_()[centre]=true;

    nr_   ()[centre] = o;
    shift_()[centre] = 1.0;
    // TODO: compute shift versus zero velocity!
    n_tot_() = (last_()+1) - first_();

    std::cout<<"before forward setup"<<std::endl;

    //Now is the perfect time to setup the boundary conditions and data for the forward ray
    comoving_ray_bdy_setup_forward<approx>(model);

    std::cout<<"dtau first+1: "<<delta_tau_()(first_()+1, 0)<<std::endl;
    //hmm, seem to be correct for now

    std::cout<<"after forward setup"<<std::endl;

    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[first_()])==rayidx)
    {
        std::cout<<"setting bdy intensities at first"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[first_()],k) += lspec.quadrature.weights[z] * wt * intensities(first_(), freqid);// Su_()[centre];
            //TODO: compute lambda term
        }
    }
    else
    {
        std::cout<<"not setting bdy intensities at first"<<std::endl;
    }

    std::cout<<"after forward boundary thing"<<std::endl;

    Size rayposidx=first_()+1;//ray position index -> point index through nr[rayposidx]
    //from first_+1 to last_ index, trace the ray in
    std::cout<<"last_(): "<<last_()<<std::endl;
    while (rayposidx<=last_())
    {
        const Real shift_next=2.0-shift_()[rayposidx];
        const Real shift_curr=2.0-shift_()[rayposidx-1];
        const bool is_upward_disc=(shift_next>=shift_curr);
        std::cout<<"solving single step"<<std::endl;
        std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        std::cout<<"dtau rayposidx: "<<delta_tau_()(rayposidx, 0)<<std::endl;
        solve_comoving_single_step (model, rayposidx, rayidx, r, is_upward_disc, true);
        rayposidx++;
    }

    std::cout<<"after forward iterations"<<std::endl;

    //Same procedure for the backward ray

    //Now is the perfect time to setup the boundary conditions and data for the backward ray
    comoving_ray_bdy_setup_backward<approx>(model);

    std::cout<<"after backwards setup"<<std::endl;

    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[last_()])==rayidx)
    {
        std::cout<<"setting bdy intensities at last"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[last_()],k) += lspec.quadrature.weights[z] * wt * intensities(last_(), freqid);// Su_()[centre];
            //TODO: compute lambda term
        }
    }


    std::cout<<"after backwards boundary thing"<<std::endl;

    rayposidx=last_()-1;//ray position index -> point index through nr[rayposidx]
    //from last_()-1 to first_ index, trace the ray backwards
    while (rayposidx>=first_())
    {
        const Real shift_next=shift_()[rayposidx];
        const Real shift_curr=shift_()[rayposidx+1];
        std::cout<<"curr ray index: "<<rayposidx+1<<std::endl;
        std::cout<<"next ray index: "<<rayposidx<<std::endl;
        const bool is_upward_disc=(shift_next>=shift_curr);
        solve_comoving_single_step (model, rayposidx, rayidx, r, is_upward_disc, false);
        rayposidx--;
    }

    std::cout<<"first ray index: "<<first_()<<std::endl;
    std::cout<<"last ray index: "<<last_()<<std::endl;
    std::cout<<"total number points: "<<n_tot_()<<std::endl;

    std::cout<<"after backwards iteration"<<std::endl;
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
//We also need the shift compared to the static frame
//Even though I disagree on the shift computation direction, other function such as get_eta_and_chi rely on this specific behavior. Ergo I must change it on a higher level.
// template <Frame frame>
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
        double shift_crt = geometry.get_shift <Rest> (o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift <Rest> (o, r, nxt, Z  );

        set_data_indicate_point (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, crt, Z, dZ);
            shift_nxt = geometry.get_shift <Rest> (o, r, nxt, Z    );

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
        //The last point of the interpolation segment is evidently a real point
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
