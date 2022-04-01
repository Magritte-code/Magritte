#include <set> //for checking whether we already have computed the intensity for a given point

template <Frame frame>
inline void Solver :: setup (Model& model)
{
    const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size  width = model.parameters.nfreqs();
    const Size  n_o_d = model.parameters.n_off_diag;

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


//FIXME: check what we now need
// //For the comoving formulation, we need to do some extra setup
// inline void Solver :: comoving_setup (Model& model)
// {
//     //compute a lower bound for dshift_max
//     Real dshift_max_lower_bound = std::numeric_limits<Real>::max();
//     Real min_inverse_mass = std::numeric_limits<Real>::max();
//
//     for (const LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
//     {
//         const Real new_inverse_mass   = lspec.linedata.inverse_mass;
//
//         if (min_inverse_mass > new_inverse_mass)
//         {
//             min_inverse_mass = new_inverse_mass;
//         }
//     }
//
//     for (Size o=0; o<model.parameters.npoints(); o++)
//     {
//         const Real new_dshift_max = model.parameters.max_width_fraction
//                                     * model.thermodynamics.profile_width (min_inverse_mass, o);
//
//         if (dshift_max_lower_bound > new_dshift_max)
//         {
//             dshift_max_lower_bound = new_dshift_max;
//         }
//     }
//
//     // for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
//     // {
//     //
//     // }
// }
template <Frame frame>
inline void Solver :: static_setup (Model& model)
{
    // const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size hnrays = model.parameters.hnrays();

    length     = 2 * get_ray_lengths_max <frame> (model) + 1;
    width = model.parameters.nfreqs();
    centre     = length/2;

    points_to_trace_ray_through.resize(hnrays);
    for (Size i = 0; i<hnrays; i++)
    {
        points_to_trace_ray_through[i].resize(model.parameters.npoints());//probably way too much, needs to be resized back down
    }
    // Matrix<Size> n_rays_through_point;
    n_rays_through_point.resize(model.parameters.hnrays(), model.parameters.npoints());
    // n_points_to_trace_ray_through.resize(model.parameters.hnrays());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        eta_c_       (i).resize (width);
        eta_n_       (i).resize (width);

        chi_c_       (i).resize (width);
        chi_n_       (i).resize (width);

        // inverse_chi_ (i).resize (length);

        tau_         (i).resize (width);

        real_pt_(i).resize(length);
        temp_intensity_(i).resize(width);
    }

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
        const Real new_dshift_max = model.parameters.max_width_fraction
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
    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        accelerated_for (o, model.parameters.npoints(),
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
    const Size      rsav)
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, rdir, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        n_rays_through_point(rsav,nxt)++;

        Size         crt = o;

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;

                  nxt = geometry.get_next(o, rdir, nxt, Z, dZ);

                  n_rays_through_point(rsav,nxt)++;
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
            if (n_rays_through_point(rr,o)>0)
            {continue;}

            //seemingly no ray has been traced through this point, so we must trace a ray through it
            points_to_trace_ray_through[rr][n_rays_to_trace]=o;
            n_rays_to_trace++;
            //tracing rays is symmetric, so only keep for the first half of the ray directions

            //trace ray through point
            n_rays_through_point(rr,o)++;
            //antipod has exactly same number of rays through the point, so do not save

            //now trace rest of rays
            trace_ray_points(model.geometry, o, rr, rr);
            trace_ray_points(model.geometry, o, ar, rr);
        }
        // n_points_to_trace_ray_through[rr]=n_rays_to_trace;
        points_to_trace_ray_through[rr].resize(n_rays_to_trace);//and now the correct size, instead of parameters.npoints()

    })

    //debug print stuff
    for (Size rr=0; rr<model.parameters.hnrays(); rr++)
    {
        std::cout<<"rr: "<<rr<<" size points_to_trace_ray_through: "<<points_to_trace_ray_through[rr].size()<<std::endl;
        for (Size idx=0; idx<points_to_trace_ray_through[rr].size(); idx++)
        {
            std::cout<<"point: "<<points_to_trace_ray_through[rr][idx]<<std::endl;
        }
        std::cout<<"number of rays per point"<<std::endl;
        for (Size p=0; p<model.parameters.npoints(); p++)
        {
            std::cout<<"point: "<<p<<"#: "<<n_rays_through_point(rr, p)<<std::endl;
        }

    }


}


// inline void Solver :: trace (Model& model)
// {
//     for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
//     {
//         const Size ar = model.geometry.rays.antipod[rr];
//
//         cout << "rr = " << rr << endl;
//
//         accelerated_for (o, model.parameters.npoints(),
//         {
//             const Real dshift_max = get_dshift_max (model, o);
//             // const Real dshift_max = 1.0e+99;
//
//             model.geometry.lengths[model.parameters.npoints()*rr+o] =
//                 trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, +1, centre+1, centre+1) + 1
//               - trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, -1, centre+1, centre  );
//         })
//
//         pc::accelerator::synchronize();
//     }
//
//     model.geometry.lengths.copy_ptr_to_vec ();
// }


inline void Solver :: solve_shortchar_order_0 (Model& model)
{
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    model.radiation.initialize_J();

    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters.npoints(),
        {
            // const Real dshift_max = get_dshift_max (o);
            const Real dshift_max = 1.0e+99;

            solve_shortchar_order_0 (model, o, rr, dshift_max);
            solve_shortchar_order_0 (model, o, ar, dshift_max);

            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}

//TODO: determine order of solver & implement multiple version (or just highest-order one)
inline void Solver :: solve_shortchar_static(Model& model)
{
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    model.radiation.initialize_J();
    model.radiation.initialize_I();//also necessary due to averaging intensities over different rays
    model.radiation.initialize_u();//also necessary due to partially filling in u

    //FIXME: also zero out all intensities
    accelerated_for (rr, model.parameters.hnrays(),
    {
    // for (Size rr = 0; rr < model.parameters.hnrays(); rr++)//parallelize over ray directions!; but first check what results we wish to write independently
    // {
        const Size ar = model.geometry.rays.antipod[rr];//antipod is probably not to useful, as we will probably be computing for both ray and antipod at the same time (as we need to trace the entire ray)

        cout << "--- rr = " << rr << endl;

        // std::vector<bool> points_already_computed(model.parameters.npoints(), false);

        //We have precomputed which points to trace through; with some smart collection of results, we might be able to parallelize this (collecting intensity at the end)
        //Eg threadprivate J, which are at the end summed (not in parallel)
        for (Size i=0; i<points_to_trace_ray_through[rr].size(); i++)
        // for (Size i=0; i<n_points_to_trace_ray_through[rr]; i++)
        {
            const Size o=points_to_trace_ray_through[rr][i];

            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;
            real_pt_()[centre]= true;

            first_() = trace_ray_static <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray_static <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                //Move dZ values, such that we have a continguous array of nonzero values (as dz[centre]=0)
                // for (Size idx=centre; idx<last_(); idx++)
                // {dZ_()[idx]=dZ_()[idx+1];}
                // Never mind, dz is apparently already a continguous array. Makes sense, as the feautrier solver also has no special exceptions

                solve_shortchar_static(model, o, rr, dshift_max);
                //FIXME: add threadprivate diagonal lambda thingy (but in the solver ofcourse)
                // then sum at the end over all results in all threads
                // for (Size f = 0; f < model.parameters.nfreqs(); f++)
                // {
                //     // solve_feautrier_order_2 (model, o, rr, ar, f);
                //
                //     model.radiation.u(rr,o,f)  = Su_()[centre];
                //     model.radiation.J(   o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];
                //
                //     update_Lambda (model, rr, f);
                // }
            }
            else
            {
                for (Size f = 0; f < model.parameters.nfreqs(); f++)
                {
                    //normally, we should divide by the number of rays through the point, though these are boundary points; so the intensity is always given by the boundary intensity
                    model.radiation.I(rr,o,f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.I(ar,o,f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    // model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    // model.radiation.J(   o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
                }
            }
            //FIXME: for all frequencies is currently not well defined, as it should denote some other value
            //Note: in later implementation, the number of frequencies might even change per ray
            // for (Size f = 0; f < model.parameters.nfreqs(); f++)
            // {
            //     model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
            //     model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            // }
        }
    // }
    })

    // compute_u_static(model);
    //and now compute J and u
    // TODO: moved to other function
    // std::cout<<"computing J and u"<<std::endl;
    for (Size o = 0; o<model.parameters.npoints(); o++)
    {
        for (Size f=0; f<model.parameters.nfreqs(); f++)
        {
            for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
            {
                // const Size ar = model.geometry.rays.antipod[rr];
                // model.radiation.u(rr,o,f)  = half * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                // model.radiation.J(   o,f) += model.geometry.rays.weight[rr] * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                model.radiation.J(   o,f) += 2.0 * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
            }
        }
    }
    // std::cout<<"finished computing"<<std::endl;

}

// For the static frame formulation, we need to shift the frequency because we need to fill in the comoving intensities
// should be computed at the end of the ray computation? No, it should be computed during computing the ray intensities
inline void Solver :: compute_u_static(Model& model)
{
    //for all certainty, put u to zero
    //first precompute vector with all frequencies
    // .resize(parameters.nfreqs())
    for (Size i=first_(); i<=last_(); i++)
    {
        Size p=nr_()[i];//point index
        Real shift=shift_()[i];//the doppler shift
        for (Size rr=0; rr<model.parameters.hnrays(); rr++)
        {
            const Size ar = model.geometry.rays.antipod[rr];
            //TODO: replace with new definition; which is slightly different for spherically symmetric grids
            // const Real rr_shift=2-model.geometry.get_shift <Rest> (p, rr, p, 0.0);
            const Real rr_shift=shift;
            // shift of other ray direction should use opposite v, so 2-(1-v)
            const Real ar_shift=2-rr_shift;
            // Size lower_freq_rayidx=0
            // Size higher_freq_rayidx=

            //first get the current shift
            //For all frequencies, shift the static intensities
            // Real temp_intensity=model.radiation.I(rr,o,f);
            // if (shift<1)//ray rr freqs should move to the left
            // {
            Real leftmost_static_freq=model.radiation.frequencies.nu(p,0);
            Real rightmost_static_freq=model.radiation.frequencies.nu(p,model.parameters.nfreqs()-1);
            Size rr_static_freq_idx=0;
            Size ar_static_freq_idx=0;

            Real rr_interpolated_I=0;
            Real ar_interpolated_I=0;
            for (Size f=0; f<model.parameters.nfreqs(); f++)
            {
                const Real rr_shifted_freq=model.radiation.frequencies.nu(p,f)*rr_shift;
                const Real ar_shifted_freq=model.radiation.frequencies.nu(p,f)*ar_shift;
                //so first for the the ray with index rr
                if (rr_shifted_freq<=leftmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,0);
                    rr_interpolated_I=model.radiation.I(rr,p,0);
                    // continue;
                }
                else if (rr_shifted_freq>=rightmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                  // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,model.parameters.nfreqs()-1);
                    rr_interpolated_I=model.radiation.I(rr,p,model.parameters.nfreqs()-1);
                    // continue;
                }
                else
                {
                    Real right_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx+1);
                    //otherwise we can just start iterating over the vector
                    while (rr_shifted_freq>right_static_freq)
                    {
                        //iterate until next right freq is larger than shifted freq
                        rr_static_freq_idx++;
                        right_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx+1);
                    }
                    Real left_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx);
                    Real left_static_intensity=model.radiation.I(rr,p,rr_static_freq_idx);
                    Real right_static_intensity=model.radiation.I(rr,p,rr_static_freq_idx+1);
                    //and now interpolate between the left and right static intensities
                    rr_interpolated_I=interpolate_linear(left_static_freq, left_static_intensity, right_static_freq, right_static_intensity, rr_shifted_freq);
                }

                //and also for the antipodal ray
                if (ar_shifted_freq<=leftmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,0);
                    ar_interpolated_I=model.radiation.I(ar,p,0);
                    // continue;
                }
                else if (ar_shifted_freq>=rightmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                  // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,model.parameters.nfreqs()-1);
                    ar_interpolated_I=model.radiation.I(ar,p,model.parameters.nfreqs()-1);
                    // continue;
                }
                else
                {
                    Real right_static_freq=model.radiation.frequencies.nu(p,ar_static_freq_idx+1);
                    //otherwise we can just start iterating over the vector
                    while (ar_shifted_freq>right_static_freq)
                    {
                        //iterate until next right freq is larger than
                        ar_static_freq_idx++;
                        right_static_freq=model.radiation.frequencies.nu(p,ar_static_freq_idx+1);
                    }
                    Real left_static_freq=model.radiation.frequencies.nu(p,ar_static_freq_idx);
                    Real left_static_intensity=model.radiation.I(rr,p,ar_static_freq_idx);
                    Real right_static_intensity=model.radiation.I(ar,p,ar_static_freq_idx+1);
                    //and now interpolate between the left and right static intensities
                    ar_interpolated_I=interpolate_linear(left_static_freq, left_static_intensity, right_static_freq, right_static_intensity, ar_shifted_freq);
                }

                model.radiation.u(rr,p,f) = half * (rr_interpolated_I+ar_interpolated_I);
                    // model.radiation.J(   o,f) += 2.0 * model.geometry.rays.weight[rr] *   model.radiation.u(rr,o,f);
                    // model.radiation.u(rr,o,f) = half * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                    // model.radiation.J(   o,f) += 2.0 * model.geometry.rays.weight[rr] *   model.radiation.u(rr,o,f);
                // }
            }
        }
    }
}


accel inline void Solver :: solve_shortchar_static (
          Model& model,
    const Size   o,
    const Size   rr,
    // const Size   ar, err, can be computed on the fly
    const double dshift_max)
{
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();
    // Vector<Real>& dZ = dZ_();
    Vector<double>& shift = shift_();

    //TODO
    // Vector<Size>& n_rays_through_point = n_rays_through_point[];//number of rays through a point
    Vector<unsigned char>& real_pt= real_pt_();
    //otherwise use char and check char equality // or other very short data type

    Vector<Real>& tau = tau_();
    Vector<Size>& nr = nr_();

    Vector<Real>& temp_intensity=temp_intensity_();
    // //zero it out before usage
    // for (Size f = 0; f<model.parameters.nfreqs(); f++)
    // {
    //     temp_intensity[f]=0.0;
    // }

    // temp_intensity.resize(model.parameters.nfreqs());


    // double  Z = 0.0;   // distance along ray
    // double dZ = 0.0;   // last distance increment
    // std::cout<<"set boundary I"<<std::endl;

    //note: boundary value can always be filled in
    for (Size f = 0; f < model.parameters.nfreqs(); f++)
    {
        model.radiation.I(rr,nr[first_()],f) = boundary_intensity(model, nr[first_()], model.radiation.frequencies.nu(o, f));
        //Shift is in the wrong direction
        model.radiation.u(rr,nr[first_()],f)+= half * boundary_intensity(model, nr[first_()], model.radiation.frequencies.nu(o, f) * (2.0-shift[first_()]))/n_rays_through_point(rr,first_()); //as we have the boundary value, why not directly fill in the comoving freq
        temp_intensity[f]=model.radiation.I(rr,nr[first_()],f);
    }

    // std::cout<<"solving for ray"<<std::endl;
    // std::cout<<"first_: "<<first_()<<std::endl;
    // std::cout<<"last_: "<<last_()<<std::endl;

    //first solve starting from direction rr
    for (Size idx=first_()+1; idx<=last_(); idx++)
    {
        const Size curr_point = nr[idx-1];
        const Size next_point = nr[idx];
        // std::cout<<"curr point: "<<curr_point<<std::endl;
        // std::cout<<"next_point: "<<next_point<<std::endl;
        // std::cout<<"nr_"<<nr[idx]<<std::endl;
        const double dZ = dZ_()[idx-1];//dZ_ contains nonzero values from first() to last()-1
        // std::cout<<"dz: "<<dZ<<std::endl;


        // Shift should be turned the other way around, so 2-(1-v)=1+v
        const double shift_curr=2.0-shift[idx-1];
        const double shift_next=2.0-shift[idx];
        // const double shift_curr=shift[idx-1];
        // const double shift_next=shift[idx];
        std::cout<<"curr point: "<<curr_point<<std::endl;
        std::cout<<"shift curr: "<<shift_curr<<std::endl;
        // std::cout<<"shift next: "<<shift_next<<std::endl;

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            get_eta_and_chi_static (model, curr_point, freq, shift_curr, eta_c[f], chi_c[f]);
            get_eta_and_chi_static (model, next_point, freq, shift_next, eta_n[f], chi_n[f]);
            // get_eta_and_chi (model, curr_point, freq*shift_curr, eta_c[f], chi_c[f]);
            // get_eta_and_chi (model, next_point, freq*shift_next, eta_n[f], chi_n[f]);


            const Real Scurr=eta_c[f]/chi_c[f];
            const Real Snext=eta_n[f]/chi_n[f];

            const Real dtau = trap (chi_c[f], chi_n[f], dZ);

            // both first order accurate
            // const Real source_term = -expm1(-dtau)*(Scurr+Snext)/2.0;//TODO: add variations
            // const Real source_term = -expm1(-dtau/2.0)*(exp(-dtau/2.0)*Scurr+Snext);//TODO: add variations

            // second order accurate
            const Real common_term = (dtau+expm1(-dtau))/dtau;
            const Real source_term = common_term*Scurr+(-common_term-expm1(-dtau))*Snext;


            // const Real onemexpmintau = -expm1(-dtau);
            // std::cout<<"chi_c: "<<chi_c[f]<<std::endl;
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"expm1(-dtau):"<<expm1(-dtau)<<std::endl;
            //
            // std::cout<<"source part: "<<onemexpmintau*source_term<<std::endl;
            temp_intensity[f]=temp_intensity[f]*exp(-dtau)+source_term;

            // model.radiation.I(rr, next_point, f)=model.radiation.I(rr, curr_point, f)*exp(-dtau)+onemexpmintau*source terms;
        }
        std::cout<<"n_rays_through_point"<<n_rays_through_point(rr,next_point)<<std::endl;
        if (real_pt[idx])
        {
            // std::cout<<"real_pt idx: "<<idx<<std::endl;
            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                //average over all rays which are traced through the point
                model.radiation.I(rr, next_point, f)+=temp_intensity[f]/n_rays_through_point(rr,next_point);
                //DEBUG: just set to intensity of first point
                // model.radiation.I(rr, next_point, f)+=temp_intensity[0]/n_rays_through_point(rr,next_point);
            }

            //while we still have computed the shift, we can fill in the u; TODO: this is inefficient, but handy for bugfixing
            const Real leftmost_static_freq=model.radiation.frequencies.nu(next_point,0);
            const Real rightmost_static_freq=model.radiation.frequencies.nu(next_point,model.parameters.nfreqs()-1);
            // std::cout<<"leftmost_static_freq: "<<leftmost_static_freq<<std::endl;
            // std::cout<<"rightmost_static_freq: "<<rightmost_static_freq<<std::endl;

            Size static_freq_idx=0;
            Real interpolated_I=0;
            for (Size f=0; f<model.parameters.nfreqs(); f++)
            {
                const Real shifted_freq=model.radiation.frequencies.nu(next_point,f)*shift_next;
                // std::cout<<"shifted freq: "<<shifted_freq<<std::endl;
                //so first for the the ray with index rr
                if (shifted_freq<=leftmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,0);
                    // std::cout<<"path: less than leftmost"<<std::endl;
                    interpolated_I=temp_intensity[0];
                    // continue;
                }
                else if (shifted_freq>=rightmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,model.parameters.nfreqs()-1);
                    // std::cout<<"path: more than rightmost"<<std::endl;
                    interpolated_I=temp_intensity[model.parameters.nfreqs()-1];
                    // continue;
                }
                else
                {
                    Real right_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx+1);
                    //otherwise we can just start iterating over the vector
                    while (shifted_freq>right_static_freq)
                    {
                        //iterate until next right freq is larger than shifted freq
                        static_freq_idx++;
                        right_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx+1);
                    }
                    Real left_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx);
                    Real left_static_intensity=temp_intensity[static_freq_idx];
                    Real right_static_intensity=temp_intensity[static_freq_idx+1];
                    //and now interpolate between the left and right static intensities
                    interpolated_I=interpolate_linear(left_static_freq, left_static_intensity, right_static_freq, right_static_intensity, shifted_freq);
                }
                model.radiation.u(rr, next_point, f)+=half * interpolated_I/n_rays_through_point(rr,next_point);
            }
                //This is somewhat more complicated, due to needing to interpolate the current intensity values
            // }
        }

        // //while we still have computed the shift, we can fill in the u; TODO: this is inefficient, but handy for bugfixing
        // Real leftmost_static_freq=model.radiation.frequencies.nu(p,0);
        // Real rightmost_static_freq=model.radiation.frequencies.nu(p,model.parameters.nfreqs()-1);
        // Size static_freq_idx=0;
        //
        // Real interpolated_I=0;
        // for (Size f=0; f<model.parameters.nfreqs(); f++)
        // {
        //     const Real rr_shifted_freq=model.radiation.frequencies.nu(p,f)*rr_shift;
        //     const Real ar_shifted_freq=model.radiation.frequencies.nu(p,f)*ar_shift;
        //     //so first for the the ray with index rr
        //     if (shifted_freq<=leftmost_static_freq)
        //     {//cant really fill in anything else except maybe the boundary value; TODO
        //         // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,0);
        //         rr_interpolated_I=model.radiation.I(rr,p,0);
        //         // continue;
        //     }
        //     else if (rr_shifted_freq>=rightmost_static_freq)
        //     {//cant really fill in anything else except maybe the boundary value; TODO
        //       // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,model.parameters.nfreqs()-1);
        //         rr_interpolated_I=model.radiation.I(rr,p,model.parameters.nfreqs()-1);
        //         // continue;
        //     }
        //     else
        //     {
        //         Real right_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx+1);
        //         //otherwise we can just start iterating over the vector
        //         while (rr_shifted_freq>right_static_freq)
        //         {
        //             //iterate until next right freq is larger than shifted freq
        //             rr_static_freq_idx++;
        //             right_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx+1);
        //         }
        //         Real left_static_freq=model.radiation.frequencies.nu(p,rr_static_freq_idx);
        //         Real left_static_intensity=model.radiation.I(rr,p,rr_static_freq_idx);
        //         Real right_static_intensity=model.radiation.I(rr,p,rr_static_freq_idx+1);
        //         //and now interpolate between the left and right static intensities
        //         rr_interpolated_I=interpolate_linear(left_static_freq, left_static_intensity, right_static_freq, right_static_intensity, rr_shifted_freq);
        //     }
        //
        //     model.radiation.u(rr,p,f) += half * interpolated_I;

    }

    // std::cout<<"Solving for reverse ray"<<std::endl;
    //Reverse direction
    const Size ar = model.geometry.rays.antipod[rr];

    //note: boundary value can always be filled in
    for (Size f = 0; f < model.parameters.nfreqs(); f++)
    {
        model.radiation.I(ar,nr[last_()],f) = boundary_intensity(model, nr[last_()], model.radiation.frequencies.nu(o, f));
        temp_intensity[f]=model.radiation.I(ar,nr[last_()],f);
        model.radiation.u(rr,nr[last_()],f)+= half * boundary_intensity(model, nr[last_()], model.radiation.frequencies.nu(o, f) * shift[last_()])/n_rays_through_point(ar,last_()); //as we have the boundary value, why not directly fill in the comoving freq
    }

    //and also solve from the reverse direction
    for (Size idx=last_(); idx>first_(); idx--)//index manipulation, as Size is unsigned
    {
        const Size curr_point = nr[idx];
        const Size next_point = nr[idx-1];
        // std::cout<<"curr point: "<<curr_point<<std::endl;
        // std::cout<<"next_point: "<<next_point<<std::endl;
        // std::cout<<"nr_"<<nr[idx]<<std::endl;
        const double dZ = dZ_()[idx-1];//dZ_ contains nonzero values from first() to last()-1
        // std::cout<<"dz: "<<dZ<<std::endl;

        // Shift should be turned the other way around in the reverse direction, so 2-(1-v)=1+v
        // const double shift_curr=2.0-shift[idx];
        // const double shift_next=2.0-shift[idx-1];
        const double shift_curr=shift[idx];
        const double shift_next=shift[idx-1];
        // std::cout<<"rev shift curr: "<<shift_curr<<std::endl;
        // std::cout<<"rev shift next: "<<shift_next<<std::endl;


        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            get_eta_and_chi_static (model, curr_point, freq, shift_curr, eta_c[f], chi_c[f]);
            get_eta_and_chi_static (model, next_point, freq, shift_next, eta_n[f], chi_n[f]);
            // get_eta_and_chi (model, curr_point, freq*shift_curr, eta_c[f], chi_c[f]);
            // get_eta_and_chi (model, next_point, freq*shift_next, eta_n[f], chi_n[f]);

            const Real Scurr=eta_c[f]/chi_c[f];
            const Real Snext=eta_n[f]/chi_n[f];

            const Real dtau = trap (chi_c[f], chi_n[f], dZ);

            // both first order accurate
            // const Real source_term = -expm1(-dtau)*(Scurr+Snext)/2.0;//TODO: add variations
            // const Real source_term = -expm1(-dtau/2.0)*(exp(-dtau/2.0)*Scurr+Snext);//TODO: add variations

            // second order accurate
            const Real common_term = (dtau+expm1(-dtau))/dtau;
            const Real source_term = common_term*Scurr+(-common_term-expm1(-dtau))*Snext;


            //TODO: check which dZ value I need, or compute it yourself using Z


            // const Real onemexpmintau = -expm1(-dtau);
            // std::cout<<"chi_c: "<<chi_c[f]<<std::endl;
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"expm1(-dtau):"<<expm1(-dtau)<<std::endl;
            //
            // std::cout<<"source part: "<<onemexpmintau*source_term<<std::endl;
            temp_intensity[f]=temp_intensity[f]*exp(-dtau)+source_term;

            // model.radiation.I(rr, next_point, f)=model.radiation.I(rr, curr_point, f)*exp(-dtau)+onemexpmintau*source terms;
        }

        if (real_pt[idx])
        {
            // std::cout<<"real_pt idx: "<<idx<<std::endl;
            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                //average over all rays which are traced through the point
                model.radiation.I(ar, next_point, f)+=temp_intensity[f]/n_rays_through_point(rr,next_point);
                //DEBUG: just set to intensity of first point
                // model.radiation.I(ar, next_point, f)+=temp_intensity[0]/n_rays_through_point(rr,next_point);
            }

            // for (Size f = 0; f < model.parameters.nfreqs(); f++)
            // {

            //while we still have computed the shift, we can fill in the u; TODO: this is inefficient, but handy for bugfixing
            Real leftmost_static_freq=model.radiation.frequencies.nu(next_point,0);
            Real rightmost_static_freq=model.radiation.frequencies.nu(next_point,model.parameters.nfreqs()-1);
            Size static_freq_idx=0;

            Real interpolated_I=0;
            for (Size f=0; f<model.parameters.nfreqs(); f++)
            {
                const Real shifted_freq=model.radiation.frequencies.nu(next_point,f)*shift_next;
                //so first for the the ray with index rr
                if (shifted_freq<=leftmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,0);
                    interpolated_I=temp_intensity[0];
                    // continue;
                }
                else if (shifted_freq>=rightmost_static_freq)
                {//cant really fill in anything else except maybe the boundary value; TODO
                    // model.radiation.u(rr,o,f)+=model.radiation.I(rr,o,model.parameters.nfreqs()-1);
                    interpolated_I=temp_intensity[model.parameters.nfreqs()-1];
                    // continue;
                }
                else
                {
                    Real right_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx+1);
                    //otherwise we can just start iterating over the vector
                    while (shifted_freq>right_static_freq)
                    {
                        //iterate until next right freq is larger than shifted freq
                        static_freq_idx++;
                        right_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx+1);
                    }
                    Real left_static_freq=model.radiation.frequencies.nu(next_point,static_freq_idx);
                    Real left_static_intensity=temp_intensity[static_freq_idx];
                    Real right_static_intensity=temp_intensity[static_freq_idx+1];
                    //and now interpolate between the left and right static intensities
                    interpolated_I=interpolate_linear(left_static_freq, left_static_intensity, right_static_freq, right_static_intensity, shifted_freq);
                }
                model.radiation.u(rr, next_point, f)+=half * interpolated_I/n_rays_through_point(rr,next_point);
            }
            //This is somewhat more complicated, due to needing to interpolate the current intensity values
            // }
        }

        //while we still have the shift, we can compute the u we can fill in


    }
}

    // Size crt = o;
    // Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    // if (model.geometry.valid_point (nxt))
    // {
    //     double shift_c = 1.0;
    //     double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);
    //
    //     for (Size f = 0; f < model.parameters.nfreqs(); f++)
    //     {
    //         const Real freq = model.radiation.frequencies.nu(o, f);
    //
    //         get_eta_and_chi (model, crt, freq,         eta_c[f], chi_c[f]);
    //         get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);
    //
    //         const Real drho = trap (eta_c[f], eta_n[f], dZ);
    //         const Real dtau = trap (chi_c[f], chi_n[f], dZ);
    //
    //         tau[f]                   = dtau;
    //         model.radiation.I(r,o,f) = drho * expf(-tau[f]);
    //     }
    //
    //     while (model.geometry.not_on_boundary (nxt))
    //     {
    //         crt     = nxt;
    //         shift_c = shift_n;
    //           eta_c =   eta_n;
    //           chi_c =   chi_n;
    //
    //         model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_n);
    //
    //         for (Size f = 0; f < model.parameters.nfreqs(); f++)
    //         {
    //             const Real freq = model.radiation.frequencies.nu(o, f);
    //
    //             get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);
    //
    //             const Real drho = trap (eta_c[f], eta_n[f], dZ);
    //             const Real dtau = trap (chi_c[f], chi_n[f], dZ);
    //
    //             tau[f]                   += dtau;
    //             model.radiation.I(r,o,f) += drho * expf(-tau[f]);
    //         }
    //     }
    //
    //     for (Size f = 0; f < model.parameters.nfreqs(); f++)
    //     {
    //         const Real freq = model.radiation.frequencies.nu(o, f);
    //
    //         model.radiation.I(r,o,f) += boundary_intensity(model, nxt, freq*shift_n) * expf(-tau[f]);
    //         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
    //     }
    // }
    //
    // else
    // {
    //     for (Size f = 0; f < model.parameters.nfreqs(); f++)
    //     {
    //         const Real freq = model.radiation.frequencies.nu(o, f);
    //
    //         model.radiation.I(r,o,f)  = boundary_intensity(model, crt, freq);
    //         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
    //     }
    // }


// inline void Solver :: solve_shortchar_order_0_comoving (Model& model)
// {
//     for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}
//
//     model.radiation.initialize_J();
//
//     accelerated_for(rr, model.parameters.hnrays(),
//     // for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
//     {
//         const Size ar = model.geometry.rays.antipod[rr];
//
//         cout << "--- rr = " << rr << endl;
//
//         //TODO: add bookkeeping about which points we already know
//
//         // accelerated_for (o, model.parameters.npoints(),
//         for (Size o=0; o<model.parameters.npoints(); o++)
//         {
//             // const Real dshift_max = get_dshift_max (o);
//             // const Real dshift_max = 1.0e+99;
//
//             // solve_shortchar_order_0 (model, o, ar, dshift_max);
//             // solve_shortchar_order_0 (model, o, rr, dshift_max);
//
// //COPY PASTED
//             const Real dshift_max = dshift_max_lower_bound;
//
//             nr_   ()[centre] = o;
//             shift_()[centre] = 1.0;
//
//             //TODO: check how the ray tracing actually works; and what thing it stores
//             //and maybe write your own variant, as the optical depths computed in this way are not used in my formula
//             //I need some weird comoving optical depth, which have some correction factor for the doppler shift
//             first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
//             last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
//             n_tot_() = (last_()+1) - first_();
//
//             //for ray and antipod, we now just compute the intensity by our comoving formula
//             //this can be done simultaneously for all frequencies
//             //the only slightly tricky part lies in the fact that we have two different opacities to use: a static one (changes greatly with doppler shifts) and a comoving one (ingores the existence of doppler shifts)
//
//             if (n_tot_() > 1)
//             {
//                 for (Size f = 0; f < model.parameters.nfreqs(); f++)
//                 {
//                     //place with something else
//                     // solve_feautrier_order_2 (model, o, rr, ar, f);
//
//                     //luckily for us, the I, u, v and J are saved comoving! But
//                     model.radiation.u(rr,o,f)  = Su_()[centre];
//                     //FIXME: only sum J at the end after synchronization; J SHOULD BE LOCAL (as we are parallelizing over the directions)
//                     model.radiation.J(   o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];
//
//                     // update_Lambda (model, rr, f); //eh, do we have a lambda operator for the shortchar?
//                 }
//             }
//             else
//             {
//                 for (Size f = 0; f < model.parameters.nfreqs(); f++)
//                 {
//                     model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
//                     model.radiation.J(   o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
//                 }
//             }
// //COPY PASTED
//
//
//             for (Size f = 0; f < model.parameters.nfreqs(); f++)
//             {
//                 model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
//                 model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
//             }
//         }//)
//
//         // pc::accelerator::synchronize();
//     })
//
//     pc::accelerator::synchronize();
//
//
//     model.radiation.I.copy_ptr_to_vec();
//     model.radiation.J.copy_ptr_to_vec();
// }


inline void Solver :: solve_feautrier_order_2 (Model& model)
{
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    model.radiation.initialize_J();

    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters.npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                for (Size f = 0; f < model.parameters.nfreqs(); f++)
                {
                    solve_feautrier_order_2 (model, o, rr, ar, f);

                    model.radiation.u(rr,o,f)  = Su_()[centre];
                    model.radiation.J(   o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];

                    update_Lambda (model, rr, f);
                }
            }
            else
            {
                for (Size f = 0; f < model.parameters.nfreqs(); f++)
                {
                    model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.J(   o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
                }
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}


inline void Solver :: image_feautrier_order_2 (Model& model, const Size rr)
{
    Image image = Image(model.geometry, rr);


    cout << "--------------------------------------" << endl;

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters.npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                image_feautrier_order_2 (model, o, rr, ar, f);

                image.I(o,f) = two*Su_()[first_()] - boundary_intensity(model, nr_()[first_()], model.radiation.frequencies.nu(o, f));
            }
        }
        else
        {
            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                image.I(o,f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
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

//For the static solvers
template <Frame frame>
accel inline Size Solver :: trace_ray_static (
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

        set_data_static (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, crt, Z, dZ);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z    );

            set_data_static (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
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


// Static solver also needs to know which points are real points in the grid, to be able to save some data
accel inline void Solver :: set_data_static (
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
///    @param[in]  p     : in dex of the cell
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
accel inline void Solver :: get_eta_and_chi (
    const Model& model,
    const Size   p,
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    // Initialize
    eta = 0.0;
    chi = 1.0e-26;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters.nlines(); l++)
    {
        const Real diff = freq - model.lines.line[l];
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);

        // cout << "prof = " << prof << "     diff = " << diff  << "     width = " << width << endl;
        // printf("prof = %le,   diff = %le,   freq = %le,   line = %le\n", prof, diff, freq, model.lines.line[l]);


        // if (isnan(eta) || isnan(chi))
        // {
        //     cout << "emmi = " << model.lines.emissivity(p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "opac = " << model.lines.opacity   (p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "widt = " << model.lines.inverse_width (p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "prof = " << prof << "   freq = " << freq << "   diff = " << diff << endl;
        // }
    }


    // cout << "eta, chi = " << eta << "  " << chi << endl;
}

///  Getter for the emissivity (eta) and the opacity (chi) in the static frame
///    @param[in]  model : reference to model object
///    @param[in]  p     : in dex of the cell
///    @param[in]  freq  : frequency (in static frame)
///    @param[in]  shift : doppler shift for line freq
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
accel inline void Solver :: get_eta_and_chi_static (
    const Model& model,
    const Size   p,
    const Real   freq,
    const double shift,
          Real&  eta,
          Real&  chi ) const
{
    // Initialize
    eta = 0.0;
    chi = 1.0e-26;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters.nlines(); l++)
    {
        const Real diff = freq - model.lines.line[l] * shift;
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l) / shift, diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);

        // cout << "prof = " << prof << "     diff = " << diff  << "     width = " << width << endl;
        // printf("prof = %le,   diff = %le,   freq = %le,   line = %le\n", prof, diff, freq, model.lines.line[l]);


        // if (isnan(eta) || isnan(chi))
        // {
        //     cout << "emmi = " << model.lines.emissivity(p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "opac = " << model.lines.opacity   (p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "widt = " << model.lines.inverse_width (p, l) << "   p = " << p << "   l = " << l << endl;
        //     cout << "prof = " << prof << "   freq = " << freq << "   diff = " << diff << endl;
        // }
    }


    // cout << "eta, chi = " << eta << "  " << chi << endl;
}




// ///  Stores the boundary intensities for later use
// accel inline void Solver :: store_far_frequencies(bool upwind, vector<Real>& boundary_intensities, vector<Real>& boundary_frequencies, )
// {
//
// }
//
// ///  Computes the boundary intensities needed for the comoving solver; also fills in the background intensity
// // Note: currently assumes the quadrature to be static; when also implementing the change in the frequency quadrature (due to temperature differences changing the line widths), revisit this
// // MAYBE TO IMPLEMENT: TODO Also can remove the outer boundary conditions from the frequency slices //but then first redefine the slices using tuples for the ranges
// //    @param[in/out]: boundary_indices
// //    @param[in/out]: intensity_map
// accel inline void Solver :: compute_boundary_frequencies_comoving(Model& model, vector<Size>& points_on_ray, vector<vector<Size>>& boundary_indices,
//                                 std::map<Real,Real>& intensity_map, vector<vector<Size>>& frequency_index_slice, Size rayindex)
// {
//     const Size tot_n_line_transitions=parameters.nlines();
//     const Size n_slices=frequency_index_slice.size()-1;
//
//     const Size total_size=(npointsonray-1)*nslices;
//
//     boundary_indices.resize(points_on_ray.size());
//
//     // vector<Real> frequencies;//FIXME: use & in function definition insteads
//     // frequencies.resize(total_size);
//     // vector<Real> intensities;
//     // intensities.resize(total_size);
//
//     //OR we just use a vector<pair<Real, Real>> as these things should be coupled together anyway.
//     // std::map<Real, Real> frequency_bdyintensity;
//     // frequency_bdyintensity.reserve(total_size);
//
//     const Size ray_origin_point=points_on_ray[0];
//
//     for (Size point_index=1; point_index<points_on_ray.size(), point_index++)//after the first point, as this can just be filled in without issues
//     {
//         const Size point_on_ray=points_on_ray[point_index];
//         const Size prev_point_on_ray=points_on_ray[points_index-1];
//         // Real dpl_shift=model.geometry.get_shift <CoMoving> (ray_origin_point, rayindex, point_on_ray, Z);//returns 1+shift
//         //check doppler shift direction
//         for (Size slice_index=0; slice_index<frequency_index_slice[point_index].size()-1; slice_index++)
//         {
//             //first check lower part of slice
//             const Real prev_lower_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index]);
//             const Real lower_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index]);//first two freqs of slice
//
//             if (lower_outer_frequency<prev_lower_outer_frequency)
//             {
//                 //set boundary value there
//                 const Real intensity_1=boundary_intensity(model, ray_origin_point, lower_outer_frequency);
//                 intensity_map[lower_outer_frequency]=intensity_1;
//                 boundary_indices.push_back(frequency_index_slice[point_index][slice_index]);
//
//
//                 const Real prev_lower_outer_frequency2=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index]+1);
//                 const Real lower_outer_frequency2=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index]+1);//first two freqs of slice
//
//                 //now also check if the second point need some boundary value (most of the time, this should be true if the previous condition is satisfied)
//                 if (lower_outer_frequency2<prev_lower_outer_frequency2)
//                 {
//                     const Real intensity_2=boundary_intensity(model, ray_origin_point, lower_outer_frequency2);
//                     intensity_map[lower_outer_frequency2]=intensity_2;
//                     boundary_indices.push_back(frequency_index_slice[point_index][slice_index]+1);
//                 }
//             }
//
//             //then check upper part of slice
//             const Real prev_upper_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index=+1]-1);
//             const Real upper_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index+1]-1);//first two freqs of slice
//
//             if (upper_outer_frequency>prev_upper_outer_frequency)
//             {
//                 //set boundary there
//                 const Real intensity_1=boundary_intensity(model, ray_origin_point, upper_outer_frequency);
//                 intensity_map[upper_outer_frequency]=intensity_1;
//                 boundary_indices.push_back(frequency_index_slice[point_index][slice_index=+1]-1);
//
//                 const Real prev_upper_outer_frequency2=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index+1]-2);
//                 const Real upper_outer_frequency2=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index+1]-2);//first two freqs of slice
//
//                 //now also check if the second point need some boundary value (most of the time, this should be true if the previous condition is satisfied)
//                 if (upper_outer_frequency2>prev_upper_outer_frequency2)
//                 {
//                     const Real intensity_2=boundary_intensity(model, ray_origin_point, upper_outer_frequency2);
//                     intensity_map[upper_outer_frequency2]=intensity_2;
//                     boundary_indices.push_back(frequency_index_slice[point_index][slice_index=+1]-2);
//                 }
//             }
//         }
//
//     //     //if positive, compute stuff on negative side
//     //     if (dpl_shift>0)
//     //     {
//     //         for (Size slice_index=0; slice_index<frequency_index_slice.size()-1; slice_index++)
//     //         {
//     //             const Real frequency_1=model.radiation.frequencies.nu(0, frequency_index_slice[slice_index]);//first two freqs of slice
//     //             const Real frequency_2=model.radiation.frequencies.nu(0, frequency_index_slice[slice_index]+1);
//     //
//     //             Real intensity_1=boundary_intensity(model, point, frequency_1*dpl_shift);
//     //             Real intensity_2=boundary_intensity(model, point, frequency_2*dpl_shift);
//     //
//     //         // frequencies[2*slice_index]=frequency_1;
//     //         // frequencies[2*slice_index+1]=frequency_2;
//     //         //
//     //         // intensities[2*slice_index]=intensity_1;
//     //         // intensities[2*slice_index+1]=intensity_2;
//     //
//     //         frequency_bdyintensity.push_back(std::make_pair())
//     //         }
//     //     }
//     //     else//if negative, compute stuff on positive side
//     //     {
//     //         for (Size slice_index=0; slice_index<frequency_index_slice.size(); slice_index++)
//     //         {
//     //             const Real frequency_1=model.radiation.frequencies.nu(o, frequency_index_slice[slice_index+1]-1);//last two freqs of slice
//     //             const Real frequency_2=model.radiation.frequencies.nu(o, frequency_index_slice[slice_index+1]-2);
//     //
//     //             Real intensity_1=boundary_intensity(model, point, frequency_1*dpl_shift);
//     //             Real intensity_2=boundary_intensity(model, point, frequency_2*dpl_shift);
//     //
//     //             frequencies[2*slice_index]=frequency_1;
//     //             frequencies[2*slice_index+1]=frequency_2;
//     //
//     //             intensities[2*slice_index]=intensity_1;
//     //             intensities[2*slice_index+1]=intensity_2;
//     //     }
//     // }
//
//     }
//
//     //TODO: sort/heapsort the vectors here // tools/heapsort
//     //maybe also check for duplicates and remove them (do point the map to the same position)
// }
//
//
// ///  Slices the frequency domain into slice which are at least max_frequency_interval away from eachother
// //   This needs to be done for every point on the ray individually
// //   @param[in] max_frequency_interval: the maximum frequency interval between frequencies in the same slice. This value should be larger than the maximal frequency difference encountered between two successive points
// //   TODO: check if we can instead use a vector for defining this adaptively
// //   @param[in/out] edge_frequency_index: for every point on the ray, stores the edges of the different frequency regions
// //   @returns the frequency indices corresponding to the first indices of the slices
// accel inline void slice_frequency_domain(Model& model, vector<Size>& points_on_ray, Real max_frequency_interval, vector<vector<Size>>& edge_frequency_index)
// {
//     // vector<vector<Size>> edge_frequency_index;
//     edge_frequency_index.resize(points_on_ray.size());
//     for (Size point_index=1; point_index<points_on_ray.size(); point_index++)
//     {
//         Size point_on_ray=points_on_ray[point_index];
//         Real deltanu=0;
//         // vector<Size> edge_frequency_index;//contains all indices such that the stored index starts at new 'block' of frequencies
//         edge_frequency_index[point_index].reserve(parameters.nlines()+1);//assume at maximum nlines, but store the first and 'last+1' index explicitly
//         edge_frequency_index[point_index].push_back(0);
//         for (Size freqid=1; freqid<parameters.nfreqs(); freqid++)
//         {
//             deltanu=model.radiation.frequencies.nu(point_on_ray, freqid)-model.radiation.frequencies.nu(point_on_ray, freqid-1);//ASSUMES the frequency discretization at all points to be the same
//             if (deltanu>max_frequency_interval)
//             {
//                 edge_frequency_index.push_back(freqid);
//             }
//         }
//         edge_frequency_index[point_index].push_back(parameters.nfreqs());
//     }
//
//     // return edge_frequency_index;
// }
//
// //adds the shifted away from intensity to the background intensity map
// accel inline void Solver :: add_intensity_to_map(std::map<Real,Real> intensity_map, Real leftfreq, Real rightfreq, Real left_intensity, Real right_intensity)
// {
//     if (leftfreq==rightfreq)//if considered frequency range has not changed, do not add anything to the map
//     {return;}
//     //FIXME: check whether minfreq!=maxfreq
//     const Real delta_freq=leftfreq-rightfreq;
//     //get all frequency values in range [minfreq, maxfreq]
//     std::map<Real,Real>::iterator itmin,itmax;
//     itmin=intensity_map.lower_bound(leftfreq);//iterator to first key with equal/higher value
//     itmax=intensity_map.upper_bound(rightfreq);//iterator to first key with higher value
//
//     for (auto it=itmin; it<itmax; it++)
//     {
//         //interpolating the intensity linearly
//         intensity_map[it]=left_intensity+(right_intensity-left_intensity)*(intensity_map[it]-leftfreq)/delta_freq;
//     }
// }
//
// // //based on the c++ docs example for lower_bound/upper_bound
// // accel inline void Solver :: remove_frequency_range_from_map(std::map<Real,Real> intensity_map, Real minfreq, Real maxfreq)
// // {
// //     std::map<Real,Real>::iterator itmin,itmax;
// //     itmin=intensity_map.lower_bound(minfreq);//iterator to first key with equal/higher value
// //     itmax=intensity_map.upper_bound(maxfreq);//iterator to first key with higher value
// //     intensity_map.erase(itlow,itup);
// // }
//
// // accel inline Real Solver :: get_stored_intensity()
// // {
// //     //If not available, just return the CMB intensity
// // }
//
//
//
//
// //TODO: add some way to flip bits of bitmap containing which points are already computed
// accel inline void Solver :: solve_comoving_order_2 (
//       Model& model,
//       const Size o, //origin point for ray
//       const Size r, //ray direction index
//       const double dshift_max) //max doppler shift for ray
//
//       Vector<Real>& eta_c = eta_c_();//current emissivity for all freqs
//       Vector<Real>& eta_n = eta_n_();//next emissivity    for all freqs
//
//       Vector<Real>& chi_c = chi_c_();//current opacity    for all freqs
//       Vector<Real>& chi_n = chi_n_();//next opacity       for all freqs
//
//       vector<Real> prev_intensity;
//       prev_intensity.resize(parameters.nfreqs());
//       //FIXME: also initialize prev_intensity with boundary intensity
//       vector<Real> curr_intensity;
//       curr_intensity.resize(parameters.nfreqs());//TODO also use fancy threadprivate stuff
//
//
//
//       //possibly also some threadlocal version of computed points bitmap? (note:vector<bool> might be faster than bitset)
//       //Vector<Real>& computed_points = computed_points_();
//
//       //optical depth not directly used in this discretization
//
//       //we do not actually care that much about the total distance along the ray
//       //But obviously it is necessary for the raytracer
//       //but maybe save the dx (however, we might just compute it on the fly; is only a difference)
//
//       double  Z = 0.0;   // distance along ray
//       double dZ = 0.0;   // last distance increment
//
//       vector<Size> points_on_ray;
//       const Size n_points_on_ray=points_on_ray.size();
//       //trace ray and antipod
//       //get all the points which lie on the current ray (spliced together)
//       //also get number points on stitched together ray
//
//       const Real FREQ_SLICE_MAX_FREQ_DIFF=3.0;//TODO: times some line width. Can be chosen adaptively for each point in principle.//FIXME: place somewhere else, maybe in function definition
//
//       vector<vector<Size>> frequency_index_slice;
//
//       vector<Size> frequency_index_slice=slice_frequency_domain(model, points_on_ray, FREQ_SLICE_MAX_FREQ_DIFF, frequency_index_slice);
//
//       vector<vector<Size>> outer_boundary_points;//for all points on the ray, stores the frequency indices of the outer boundary points
//       outer_boundary_points.reserve(n_points_on_ray);
//       for (Size point_index=1; point_index<n_points_on_ray; point_index++)
//       {
//           outer_boundary_points[point_index].reserve(4*frequency_index_slice[point_index]);//true upper bound, as every slice can only have 2 boundary points on each side
//       }
//       std::map<Real,Real> boundary_intensity;//maps the frequency to the intensity for all thrown away values of the intensity
//
//       compute_boundary_frequencies_comoving(model, points_on_ray, outer_boundary_points, boundary_intensity, frequency_index_slice, r);
//
//
//       //first two slice points are boundaries is backward discretization
//       //inner boundary only possible if from forward to downward
//       //resulting domains should only contain successive points with same discretization without the boundary points (as they are implied)
//
//       // Vector<Real> boundary_frequencies;//the frequencies at which the boundary conditions will need to be evaluated
//       // //.resize(npoints*total number transitions)
//       // Vector<Real> boundary_intensities;//the corresponding boundary intensities; start by filling in the default CMB intensity at the edge
//       // //.resize(npoints*total number transitions)
//       // //Also define map from the transition and point to the index of this
//
//       //now sort the frequencies (and then fill in the boundary intensities)
//       //also create a map (which should be co-sorted) to find out which we actually need in the vector
//
//
//       //do add bitmap denoting whether we should store the computed intensity
//
//       //somewhere also set bits of other very large bitmap O(npoints), denoting that we have/will have computed the intensities of the points
//
//       //first setup all the boundary intensities necessary (use ray length (including ghost points))
//
//       //now compute the intensities for all frequencies, for all points on the ray
//       for (Size point_index=1; point_index<n_points_on_ray; point_index++)
//       {
//           //get doppler shift ; we actually just need the frequency diff, so no longer necessary
//           // double shift_c = 1.0;
//           // double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);//returns 1+shift
//           // double relative_shift=shift_n-shift_c;//finally just the shift
//           const Size prev_point_on_ray=points_on_ray[point_index-1];
//           const Size curr_point_on_ray=points_on_ray[point_index];
//           const Size n_slices=frequency_index_slice[point_index].size()-1;
//
//           TODO compute deltax
//
//
//           // vector<Real> temp_intensity;
//           // temp_intensity.resize(parameters.nfreqs());
//           // vector<Real> prev_intensity;//TODO: also resize and fill this in with boundary intensities
//
//
//           //for all frequencies, do the rhs step
//           //if we assume the lines to be independent (far enough away from eachother), we can split this up into the different quadratures
//           //DO FIGURE OUT: how to simply find which frequencies should be grouped together; maybe compute the max doppler shift in the entire domain and base some frequency distance off that (+ some extra to be safe)
//           //if then the line frequencies lie too close together, we should somehow indicate to our algorithm that we should have a frequency quadrature together. (but then we first need to sort the line frequencies... hmm, this will take some time to work out)
//           // for (all line species)
//           // for (all radiative transitions)//<-up until (including) this point, we should treat them independently
//           // for (all quadratures)
//
//
//           // vector<vector<bool>> using_forward_discretization;//for a given point; it stores for every slice, for every frequency in the slice, whether a foward discretization is used
//           vector<std::pair<Size, Size>> inner_ranges;//for every seperate region to solve, denotes the range of the frequency region
//           vector<bool> using_forward_discretization;//for a given range, denotes the discretization to use
//           vector<std::pair<Size,Size>> inner_boundary_tuples;//contains all inner boundary pairs
//           // using_forward_discretization.resize();
//           // allocate some estimate
//           inner_ranges.reserve(2*n_slices);//assuming no merging of line profiles, this should be enough
//           using_forward_discretization.reserve(2*n_slices);//reserve the same
//           inner_boundary_tuples.reserve(n_slices);//assuming no merging of line profiles, doppler shift together with line width change should result in only a single inflection point
//
//
//
//
//           //hmm, can't we do this preprocessing step for all frequencies simultaneously?
//           for (Size slice_index=0; slice_index<frequency_index_slice[point_index].size()-1; slice_index++) //every slice can be treated independently
//           {
//               const Size freq_slice_start=frequency_index_slice[point_index][slice_index];
//               const Size freq_slice_end=frequency_index_slice[point_index][slice_index];//actually one behind the last point of the slice (eg max parameters.nfreqs())
//               for (Size f=freq_slice_start; f<freq_slice_end; f++)
//           // for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
//           // Size f=lineProducingSpecies[ls].nr_line[point][radiative transition][quadrature index]//FIXME: check whether the corresponding frequencies are ordered in a given radiative transition
//           // for(Size f = 0; f < model.parameters.nfreqs(); f++)
//               {
//                   const Real freqprev = model.radiation.frequencies.nu(prev_point_on_ray, f);//FIXME: figure out whether this is comoving freq or not
//                   const Real freqcurr = model.radiation.frequencies.nu(curr_point_on_ray, f);//TODO
//                   //Huh, do I even care about this roundabout way?; can't we just compute how far from the line center we are (first get inverse width) and multiply the line opacity with it?
//                   //No, a priori, I do not know whether dust (or some other sources beside the current line opacity) exist in the region.
//
//                   get_eta_and_chi (model, crt, freqprev, eta_c[f], chi_c[f]);
//                   get_eta_and_chi (model, nxt, freqcurr, eta_n[f], chi_n[f]);
//
//
//
//           }
//           }
//           //end for quads
//           //end for transitions
//           //end for line species
//
//           // for (Size slice_index=0; slice_index<frequency_index_slice.size()-1; slice_index++) //every slice can be treated independently
//           // for (all freqs in slice)
//           //check which region is should have upward/downward discretization
//           for (Size slice_index=0; slice_index<frequency_index_slice[point_index].size()-1; slice_index++) //every slice can be treated independently
//           {
//               Size freq_slice_start=frequency_index_slice[point_index][slice_index];
//               Size freq_slice_end=frequency_index_slice[point_index][slice_index];//actually one behind the last point of the slice (eg max parameters.nfreqs())
//               //get first freq diff, if negative, we have an outer boundary condition on the left side; repeat for the next point as we can have up to two boundary conditions on each side
//               Real prev_lower_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index]);
//               Real lower_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index]);//first freq of slice
//               if (lower_outer_frequency<prev_lower_outer_frequency)
//               {//boundary condition, so start one later
//                   freq_slice_start+=1;
//                   prev_lower_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index]);
//                   lower_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index]);//second freq of slice
//
//                   if (lower_outer_frequency<prev_lower_outer_frequency)
//                   {//boundary condition, so start one later
//                       freq_slice_start+=1;
//                   }
//               }
//
//               Real prev_upper_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index=+1]-1);
//               Real upper_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index+1]-1);//last freq of slice
//
//               if (upper_outer_frequency>prev_upper_outer_frequency)
//               {
//                   freq_slice_end-=1;
//                   Real prev_upper_outer_frequency=model.radiation.frequencies.nu(prev_point_on_ray, frequency_index_slice[point_index][slice_index=+1]-2);
//                   Real upper_outer_frequency=model.radiation.frequencies.nu(point_on_ray, frequency_index_slice[point_index][slice_index+1]-2);//second last freq of slice
//
//                   if (upper_outer_frequency>prev_upper_outer_frequency)
//                   {
//                       freq_slice_end-=1;
//                   }
//               }
//
//               //now start with assigning continguous regions
//               Size range_start=freq_slice_start;
//               // Size range_end=freq_slice_end;
//               prev_nu=model.radiation.frequencies.nu(prev_point_on_ray, range_start);
//               curr_nu=model.radiation.frequencies.nu(point_on_ray, range_start);
//               bool curr_is_upward_discretization=(curr_nu>prev_nu);
//
//               for (Size f=freq_slice_start+1; f<freq_slice_end; f++)
//               {
//                   prev_nu=model.radiation.frequencies.nu(prev_point_on_ray, f);
//                   curr_nu=model.radiation.frequencies.nu(point_on_ray, f);
//                   bool temp_is_upward_discretization=(curr_nu>prev_nu);
//                   if (curr_is_upward_discretization!=temp_is_upward_discretization)
//                   {
//                       //TODO: check for ranges of size zero? Does this actually matter?
//                       //we need to end our current region
//                       if (curr_is_upward_discretization)//going from upward to downward means inner boundary condition
//                       {
//                           //set inner bound at f-1,f
//                           inner_ranges.push_back(std::make_pair(range_start,f-1));
//                           using_forward_discretization.push_back(true);
//                           range_start=f+1;
//                           inner_boundary_tuples.push_back(std::make_pair(f-1,f));
//                           f+=1;//skip over the second inner boundary point
//                       }
//                       else//just two seperate regions
//                       {
//                           inner_ranges.push_back(std::make_pair(range_start,f));
//                           range_start=f;
//                           using_forward_discretization.push_back(false);
//                       }
//                       curr_is_upward_discretization=temp_is_upward_discretization;//discretization direction has changed
//                   }
//               }
//               //and do not forget to put a domain at the end
//               inner_ranges.push_back(std::make_pair(range_start,freq_slice_end));
//               using_forward_discretization.push_back(curr_is_upward_discretization);
//           }
//           //end for all slices
//
//           //now compute all boundary conditions
//           //first outer boundary conditions
//           for (Size boundary_point: outer_boundary_points[point_index])
//           {
//             const Real freq=model.radiation.frequencies.nu(curr_point_on_ray, f);
//             curr_intensity[boundary_point]=boundary_intensity[freq];
//           }
//           //now compute inner boundary point pairs
//           //maybe TODO: can be vectorized, as the boundary pairs should stand right next to eachother
//           for (auto pairs: inner_boundary_tuples)
//           {
//             Size firstid=inner_boundary_tuples.first;
//             Size secondid=inner_boundary_tuples.second;
//
//             const Real prevfirstfreq=model.radiation.frequencies.nu(prev_point_on_ray, firstid);
//             const Real currfirstfreq=model.radiation.frequencies.nu(curr_point_on_ray, firstid);
//             const Real prevsecondfreq=model.radiation.frequencies.nu(prev_point_on_ray, secondid);
//             const Real currsecondfreq=model.radiation.frequencies.nu(curr_point_on_ray, secondid);
//
//             const Real delta_freq=prev_second_intensity-prev_first_intensity;//freq width of interval
//             const Real first_freq_rel_diff=(currfirstfreq-prevfirstfreq)/delta_freq;//distances for interpolation
//             const Real second_freq_rel_diff=(currsecondfreq-prevfirstfreq)/delta_freq;
//
//             const Real prev_first_intensity=prev_intensity[firstid];
//             const Real prev_second_intensity=prev_intensity[secondid];
//             const Real delta_intensity=prev_second_intensity-prev_first_intensity;
//
//             //interpolate linearly the prev_intensity to the next freq values
//             Real temp_first_intensity=prev_first_intensity+delta_intensity*first_freq_rel_diff;
//             Real temp_second_intensity=prev_first_intensity+delta_intensity*second_freq_rel_diff;
//
//             //and now apply static solver to the interpolated values
//             curr_intensity[firstid]=(temp_first_intensity*(1.0-deltax*chi_c[firstid]/2.0)+deltax*(eta_c[firstid]+eta_n[firstid])/2.0)
//                                     /(1+deltax*chi_c[firstid]/2.0);
//             curr_intensity[secondid]=(temp_second_intensity*(1.0-deltax*chi_c[secondid]/2.0)+deltax*(eta_c[secondid]+eta_n[secondid])/2.0)
//                                      /(1+deltax*chi_c[secondid]/2.0);
//           }
//
//           //for all regions, we can finally start computing stuff
//           //both explicit and implicit parts can now be computed
//           for (Size range_id=0; range_id<inner_ranges.size(); range_id++)
//           {
//               Size range_start=inner_ranges[range_id].first;
//               Size range_end=inner_ranges[range_id].end;//stores index+1s
//               bool is_forward_discretization=using_forward_discretization[range_id];
//
//               if (is_forward_discretization)
//               {
//                   //we know the boundary conditions on the top side, so start iteration from there
//                   for(Size i=range_start; i<range_end; i++)//some index manipulation necessary to iterate starting from the end of the range
//                   {
//                       const Size f=range_end-1-i+range_start;
//
//                       const Real nu_0=model.radiation.frequencies.nu(prev_point_on_ray, f);
//                       const Real nu_1=model.radiation.frequencies.nu(prev_point_on_ray, f+1);
//                       const Real nu_2=model.radiation.frequencies.nu(prev_point_on_ray, f+2);
//
//                       const Real nu_diff=model.radiation.frequencies.nu(curr_point_on_ray, f)-nu_0;
//                       //compute constants for the frequency derivative
//                       //nu1 should be the next freq, nu2 should be the second next freq
//                       const Size freq_coef_2=(nu1-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu2-nu0,2));
//                       const Size freq_coef_1=(nu2-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu1-nu0,2));
//                       const Size freq_coef_0=-freq_coef_1-freq_coef_2;
//
//                       //explicit part
//                       curr_intensity[f]=prev_intensity[f]*(1.0-deltax*chi_c[f]/2.0)+deltax*(eta_c[f]+eta_n[f])/2.0
//                                         +deltax/2.0*nu_diff*(freq_coef_0*prev_intensity[f]+freq_coef_1*prev_intensity[f+1] +freq_coef_2*prev_intensity[f+2]);//moving the frequencies in the positive direction, so using the next frequencies on the positive side for estimating the derivative
//
//                       //implicit part
//                       curr_intensity[f]=(curr_intensity[f]-deltax*nu_diff/2.0*(freq_coef_1*curr_intensity[f+1]+freq_coef_2*curr_intensity[f+2]))
//                                         /(1.0+chi_n*deltax/2.0-deltax*relative_shift/2.0*freq_coef_0);
//                   }
//                   const Real prev_leftmost_freq=model.radiation.frequencies.nu(prev_point_on_ray, range_start);
//                   const Real curr_leftmost_freq=model.radiation.frequencies.nu(curr_point_on_ray, range_start);
//
//                   const Real prev_leftmost_intensity=prev_intensity[range_start];
//                   const Real curr_leftmost_intensity=curr_intensity[range_start];
//                   //now do bookkeeping for the background intensity
//                   //this region has shifted away to higher frequencies, so we need to replace the left intensities
//                   add_intensity_to_map(boundary_intensity, prev_leftmost_freq, curr_leftmost_freq, prev_leftmost_intensity, curr_leftmost_intensity);
//               }
//               else
//               {
//                   //we know the boundary conditions on the bottom side, so start iteration from there
//                   for(Size f=range_start; f<range_end; f++)
//                   {
//                       const Real nu_0=model.radiation.frequencies.nu(prev_point_on_ray, f);
//                       const Real nu_1=model.radiation.frequencies.nu(prev_point_on_ray, f-1);
//                       const Real nu_2=model.radiation.frequencies.nu(prev_point_on_ray, f-2);
//
//                       const Real nu_diff=model.radiation.frequencies.nu(curr_point_on_ray, f)-nu_0;
//
//                       //compute constants for the frequency derivative
//                       //nu1 should be the previous freq, nu2 should be the second previous freq
//                       const Size freq_coef_2=(nu1-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu2-nu0,2));
//                       const Size freq_coef_1=(nu2-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu1-nu0,2));
//                       const Size freq_coef_0=-freq_coef_1-freq_coef_2;
//
//                       curr_intensity[f]=prev_intensity[f]*(1.0-deltax*chi_c[f]/2.0)+deltax*(eta_c[f]+eta_n[f])/2.0
//                                   +deltax*nu_diff*(freq_coef_0*prev_intensity[f]+freq_coef_1*prev_intensity[f-1] +freq_coef_2*prev_intensity[f-2]);//moving the frequencies in the positive direction, so using the next frequencies on the positive side for estimating the derivative
//
//                       curr_intensity[f]=(curr_intensity[f]-deltax*nu_diff/2.0*(freq_coef_1*curr_intensity[f-1]+freq_coef_2*curr_intensity[f-2]))
//                                         /(1.0+chi_n*deltax/2.0-deltax*relative_shift/2.0*freq_coef_0);
//                   }
//                   const Real prev_rightmost_freq=model.radiation.frequencies.nu(prev_point_on_ray, range_end-1);
//                   const Real curr_rightmost_freq=model.radiation.frequencies.nu(curr_point_on_ray, range_end-1);
//
//                   const Real prev_rightmost_intensity=prev_intensity[range_end-1];
//                   const Real curr_rightmost_intensity=curr_intensity[range_end-1];
//                   //now do bookkeeping for the background intensity
//                   //this region has shifted away to lower frequencies, so we need to replace the right intensities
//                   add_intensity_to_map(boundary_intensity, curr_rightmost_freq, prev_rightmost_freq, curr_rightmost_intensity, prev_rightmost_intensity);
//               }
//           }
//
//           //bookkeeping should now be done
//
//
//           // //add freq derivative stuff // also do not forget to slice the slices into the different discretization domain (to which one can stably apply the discretization scheme)
//           // if (relative_shift>0)
//           // {
//           //
//           // for (all line species)
//           // for (all radiative transitions)
//           // //solve per frequency quadrature individually, rhs stuff should be vectorizable without problem
//           // Size starting_quad_index=TODO;//FIXME: use this to offset the index for the current intensity
//           // //OR just use some multidimensional array
//           //
//           // for (all quadratures)//minus the boundary quadratures, as they get set to the boundary intensity
//           // //DO iterate starting from the boundary!
//           //
//           // //compute constants for the frequency derivative
//           // //nu1 should be the next freq, nu2 should be the second next freq
//           // const Size freq_coef_2=(nu1-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu2-nu0,2));
//           // const Size freq_coef_1=(nu2-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu1-nu0,2));
//           // const Size freq_coef_0=-freq_coef_a-freq_coef_b;
//           //
//           // curr_intensity[quad]=prev_intensity[quad]*(1.0-deltax*chi_c/2.0)+deltax*(eta_c+eta_n)/2.0
//           //             +deltax/2.0*relative_shift*(freq_coef_0*prev_intensity[quad]+freq_coef_1*prev_intensity[quad+1] +freq_coef_2*prev_intensity[quad+2]);//moving the frequencies in the positive direction, so using the next frequencies on the positive side for estimating the derivative
//           //
//           // //end for quadratures
//           //
//           // //now apply boundary conditions (boundary intensities)
//           //
//           // //FIXME: just define function/map that gets index of boundary intensity (wait, cant we just sort the intensities and also keep the shuffled indices?)
//           // //ALSO define an array with all the frequencies required for the boundary intensities
//           // curr_intensity[nquads-2]=bdy_intensity[TODOfindindex]
//           // curr_intensity[nquads-1]=bdy_intensity[TODOfindindex]
//           //
//           // for (all quadratures)
//           // //setup eigen matrix for left hand side/but also try out manually in another attempt
//           // //err, eigen does not support sparse representations of banded matrices... (as far as I know as time of writing), so backward substitution is the way to go
//           // curr_intensity[quad]=(currintensity[quad]
//           //                       -deltax*relative_shift/2.0*(freq_coef_1*prev_intensity[quad+1]+freq_coef_2*prev_intensity[quad+2]))
//           //                      /(1.0+chi_n*deltax/2.0-deltax*relative_shift/2.0*freq_coef_0);
//           //
//           // //end for quadratures
//           // //and solve the matrix/vector equation
//           //
//           // //bookkeeping: old intensities when shifted away from, should be interpolated and put in the boundary intensities (if going over a boundary intensity)
//           // //check range which the frequencies stepped over
//           // //check which indices in the sorted boundary frequency vector this corresponds to
//           // //interpolate between the intensity at the previous and the current far end linearly I(prev, nuprev)+(nubound-nuprev)*I(curr, nucurr)/(nucurr-nuprev)
//           // //WARNING: only do this if we actually need to (divisions by zero should be avoided)
//           //
//           // //end for transitions
//           // //end for line species
//           // }
//           // else
//           // {//TODO: similar for the other case relative_shift<=0
//           //
//           //
//           // for (all line species)
//           // for (all radiative transitions)
//           // //solve per frequency quadrature individually, rhs stuff should be vectorizable without problem
//           // for (all quadratures)//minus the boundary quadratures, as they get set to the boundary intensity
//           // //DO iterate starting from the boundary! (2->nquads-1)
//           //
//           // //compute constants for the frequency derivative
//           // //nu1 should be the previous freq, nu2 should be the second previous freq
//           // const Size freq_coef_2=(nu1-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu2-nu0,2));
//           // const Size freq_coef_1=(nu2-nu0)/((nu2-nu0)*(nu1-nu0)-std::pow(nu1-nu0,2));
//           // const Size freq_coef_0=-freq_coef_a-freq_coef_b;
//           //
//           // curr_intensity[quad]=prev_intensity[quad]*(1.0-deltax*chi_c)+deltax*(eta_c+eta_n)/2.0
//           //             +deltax*relative_shift*(freq_coef_0*prev_intensity[quad]+freq_coef_1*prev_intensity[quad-1] +freq_coef_2*prev_intensity[quad-2]);//moving the frequencies in the positive direction, so using the next frequencies on the positive side for estimating the derivative
//           //             //?minus sign? because using the backward discretization (coefficients should then be *-1^(nth derivative) HOWEVER, the computed coefficients already changed sign, so no issues there
//           //
//           // //explicitly split the for loops
//           //
//           // //end for quadratures
//           //
//           // //now apply boundary conditions (boundary intensities)
//           // //FIXME: just define function/map that gets index of boundary intensity (wait, cant we just sort the intensities and also keep the shuffled indices?)
//           // //ALSO define an array with all the frequencies required for the boundary intensities
//           // curr_intensity[0]=bdy_intensity[TODOfindindex]
//           // curr_intensity[1]=bdy_intensity[TODOfindindex]
//           //
//           // for (all quadratures)
//           // //setup eigen matrix for left hand side/but also try out manually in another attempt//again, eigen does seemingly not simply support banded matrices, so implementing forward substitution
//           // curr_intensity[quad]=(currintensity[quad]
//           //                       -deltax*relative_shift/2.0*(freq_coef_1*prev_intensity[quad-1]+freq_coef_2*prev_intensity[quad-2]))
//           //                      /(1.0+chi_n*deltax/2.0-deltax*relative_shift/2.0*freq_coef_0);
//           // //end for quadratures
//           // //and solve the matrix/vector equation
//           //
//           // //end for quadratures
//           // //bookkeeping: old intensities when shifted away from, should be interpolated and put in the boundary intensities (if going over a boundary intensity)
//           // //check range which the frequencies stepped over
//           // //check which indices in the sorted boundary frequency vector this corresponds to
//           // //interpolate between the intensity at the previous and the current far end linearly I(prev, nuprev)+(nubound-nuprev)*I(curr, nucurr)/(nucurr-nuprev)
//           // //WARNING: only do this if we actually need to (divisions by zero should be avoided)
//           //
//           // //end for transitions
//           // //end for line species
//           // }
//
//
//
//           //now we have new intensities
//
//
//
//
//
//           //check whether we need to store the intensities (for all frequencies at the same time)
//           if (TODOstoreintensity)
//           {
//             //store intensity for computing J
//             //Later: also setup diagonal lambda elements (partially, as they should be summed over all directions)
//           }
//
//           //EXTRA: maybe also setup the diagonal lambda elements (but first check the computation)
//           // NOTE: lspec.lambda.add_element(nr[centre], k, nr[n], L); just adds results if same index is used twice (quite useful, if we were not computing for all comoving frequencies at once?... wait, it might simplify everything a lot)
//           // or just use som std::sum -like function
//       }
//
//
//
//
//
//
// )




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

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            get_eta_and_chi (model, crt, freq,         eta_c[f], chi_c[f]);
            get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);

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

            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                const Real freq = model.radiation.frequencies.nu(o, f);

                get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);

                const Real drho = trap (eta_c[f], eta_n[f], dZ);
                const Real dtau = trap (chi_c[f], chi_n[f], dZ);

                tau[f]                   += dtau;
                model.radiation.I(r,o,f) += drho * expf(-tau[f]);
            }
        }

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f) += boundary_intensity(model, nxt, freq*shift_n) * expf(-tau[f]);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
        }
    }

    else
    {
        for (Size f = 0; f < model.parameters.nfreqs(); f++)
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
accel inline void Solver :: solve_feautrier_order_2 (
          Model& model,
    const Size   o,
    const Size   rr,
    const Size   ar,
    const Size   f  )
{
    const Real freq = model.radiation.frequencies.nu(o, f);

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
    get_eta_and_chi (model, nr[first  ], freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi (model, nr[first+1], freq*shift[first+1], eta_n, chi_n);

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
        get_eta_and_chi (model, nr[n+1], freq*shift[n+1], eta_n, chi_n);

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
accel inline void Solver :: image_feautrier_order_2 (
          Model& model,
    const Size   o,
    const Size   rr,
    const Size   ar,
    const Size   f  )
{
    const Real freq = model.radiation.frequencies.nu(o, f);

    // cout << "o = " << o << "  f = " << f << "  " << CC*(freq / model.lines.line[0] - 1.0) << endl;

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
    get_eta_and_chi (model, nr[first  ], freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi (model, nr[first+1], freq*shift[first+1], eta_n, chi_n);

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

    // cout << "FF[first] = " << FF[first] << "  dtau_n = " << dtau_n << "   chi_n = " << chi_n << "   chi_c = " << chi_c <<  endl;

    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        // Get new radiative properties
        get_eta_and_chi (model, nr[n+1], freq*shift[n+1], eta_n, chi_n);

        inverse_chi[n+1] = 1.0 / chi_n;

        term_n = eta_n * inverse_chi[n+1];
        dtau_n = half * (chi_c + chi_n) * dZ[n];

        // cout << "dtau = " << dtau_n << endl;

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        // cout << "inverse_C[" << n << "] = " << inverse_C[n] << endl;

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

    // cout << "Bl = " << Bl << "   FF[last-1] = " << FF[last-1] << "   Bl_min_Al = " << Bl_min_Al << endl;

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    for (long n = last-1; n > first; n--) // use long in reverse loops!
    {
        Su[n] += Su[n+1] * FI[n];
    }

    Su[first] += Su[first+1] * FI[first];
}





// const Real alpha  = 1.0;
// const Real alpha2 = alpha * alpha;
// const Real h_smt  = 1.0;
// const Real h_smt2 = h_smt * h_smt;
// const Real inverse_h_smt2 = 1.0 / h_smt2;
// const Real inverse_h_smt4 = inverse_h_smt2 * inverse_h_smt2;
// const Real minus_half_inverse_h_smt2 = -0.5 * inverse_h_smt2;
//
//
// accel inline Real Solver :: kernel (const Vector3D d) const
// {
//     return alpha2 * exp(minus_half_inverse_h_smt2 * d.squaredNorm());
// }
//
//
// accel inline Real Solver :: kernel (
//     const Model& model,
//     const Size   r,
//     const Size   p1,
//     const Size   p2 ) const
// {
//     const Vector3D x1 = model.geometry.points.position[p1];
//     const Vector3D x2 = model.geometry.points.position[p2];
//
//     return kernel(x1-x2);
// }
//
//
// accel inline Real Solver :: L1_kernel (
//     const Model& model,
//     const Size   r,
//     const Size   p1,
//     const Size   p2 ) const
// {
//     const Vector3D d =   model.geometry.points.position[p1]
//                        - model.geometry.points.position[p2];
//
//     const Real g = d.dot(model.geometry.rays.direction[r]) * inverse_h_smt2;
//
//     return (chi[p1] - g) * kernel(d);
// }
//
//
// accel inline Real Solver :: L2_kernel (
//     const Model& model,
//     const Size   r,
//     const Size   p1,
//     const Size   p2 ) const
// {
//     const Vector3D d =   model.geometry.points.position[p1]
//                        - model.geometry.points.position[p2];
//
//     const Real g = d.dot(model.geometry.rays.direction[r]) * inverse_h_smt2;
//
//     return (chi[p2] + g) * kernel(d);
// }
//
//
// accel inline Real Solver :: L12_kernel (
//     const Model& model,
//     const Size   r,
//     const Size   p1,
//     const Size   p2 ) const
// {
//     const Vector3D d =   model.geometry.points.position[p1]
//                        - model.geometry.points.position[p2];
//
//     const Real g = d.dot(model.geometry.rays.direction[r]) * inverse_h_smt2;
//
//     return ((chi[p1] + g)*(chi[p2] - g) + inverse_h_smt2) * kernel(d);
// }
//
//
// accel inline void Solver :: solve_kernel_method (
//           Model& model,
//     const Size   r,
//     const Size   f )
// {
//     const Real freq = model.radiation.frequencies.nu(0, f);
//
//     // Get emissivity and opacity
//     eta.resize (model.parameters.npoints());
//     chi.resize (model.parameters.npoints());
//
//     for (Size p = 0; p < model.parameters.npoints(); p++)
//     {
//         get_eta_and_chi (model, p, freq, eta[p], chi[p]);
//     }
//
//
//     // Triplets for (sparse) covariance matrix
//     vector<Triplet<Real, Size>> triplets;
//
//     VectorXr y (model.parameters.npoints() + model.parameters.nboundary());
//     VectorXr w (model.parameters.npoints() + model.parameters.nboundary());
//
//
//     for (Size b1 = 0; b1 < model.parameters.nboundary(); b1++)
//     {
//         const Size p1 = model.geometry.boundary.boundary2point[b1];
//
//         for (Size b2 = 0; b2 < model.parameters.nboundary(); b2++)
//         {
//             const Size p2 = model.geometry.boundary.boundary2point[b2];
//
//             triplets.push_back (Triplet<Real, Size> (
//                 b1,
//                 b2,
//                 kernel (model, r, p1, p2)
//             ));
//         }
//
//         for (Size p2 = 0; p2 < model.parameters.npoints(); p2++)
//         {
//             triplets.push_back (Triplet<Real, Size> (
//                 b1,
//                 p2 + model.parameters.nboundary(),
//                 L2_kernel (model, r, p1, p2)
//             ));
//         }
//
//         y[b1] = boundary_intensity (model, p1, freq);
//     }
//
//
//     for (Size p1 = 0; p1 < model.parameters.npoints(); p1++)
//     {
//         const Size i1 = p1 + model.parameters.nboundary();
//
//         for (Size b2 = 0; b2 < model.parameters.nboundary(); b2++)
//         {
//             const Size p2 = model.geometry.boundary.boundary2point[b2];
//
//             triplets.push_back (Triplet<Real, Size> (
//                 i1,
//                 b2,
//                 L1_kernel (model, r, p1, p2)
//             ));
//         }
//
//         for (Size p2 = 0; p2 < model.parameters.npoints(); p2++)
//         {
//             triplets.push_back (Triplet<Real, Size> (
//                 i1,
//                 p2 + model.parameters.nboundary(),
//                 L12_kernel (model, r, p1, p2)
//             ));
//         }
//
//         y[i1] = eta[p1];
//     }
//
//
//     SparseMatrix<Real> covariance;
//     covariance.setFromTriplets (triplets.begin(), triplets.end());
//
//     SparseLU <SparseMatrix<Real>, COLAMDOrdering<int>> solver;
//
//     cout << "Analyzing covariance matrix..." << endl;
//     solver.analyzePattern (covariance);
//     cout << "Factoring covariance matrix..." << endl;
//     solver.factorize      (covariance);
//
//     if (solver.info() != Eigen::Success)
//     {
//         throw std::runtime_error (solver.lastErrorMessage());
//     }
//
//     cout << "Inverting covariance matrix..." << endl;
//     w = solver.solve (y);
//
//
//     for (Size p1 = 0; p1 < model.parameters.npoints(); p1++)
//     {
//         model.radiation.I(r, p1, f) = 0.0;
//
//         for (Size b2 = 0; b2 < model.parameters.nboundary(); b2++)
//         {
//             const Size p2 = model.geometry.boundary.boundary2point[b2];
//
//             model.radiation.I(r, p1, f) += kernel(model, r, p1, p2) * w[b2];
//         }
//
//         for (Size p2 = 0; p2 < model.parameters.npoints(); p2++)
//         {
//             const Size i2 = p2 + model.parameters.nboundary();
//
//             model.radiation.I(r, p1, f) += L2_kernel(model, r, p1, p2) * w[i2];
//         }
//     }
//
//     return;
// }

accel inline void Solver :: set_eta_and_chi (Model& model) const
{
    model.eta.resize (model.parameters.npoints(), model.parameters.nfreqs());
    model.chi.resize (model.parameters.npoints(), model.parameters.nfreqs());

    for (Size p = 0; p < model.parameters.npoints(); p++)
    {
        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(0, f);

            get_eta_and_chi (model, p, freq, model.eta(p,f), model.chi(p,f));
        }
    }
}


accel inline void Solver :: set_boundary_condition (Model& model) const
{
    model.boundary_condition.resize (model.parameters.nboundary(), model.parameters.nfreqs());

    for (Size b = 0; b < model.parameters.nboundary(); b++)
    {
        const Size p = model.geometry.boundary.boundary2point[b];

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(0, f);

            model.boundary_condition(b,f) = boundary_intensity (model, p, freq);
        }
    }
}
