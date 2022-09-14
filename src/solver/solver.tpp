template <Frame frame>
inline void Solver :: setup (Model& model)
{
    const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size  width = model.parameters->nfreqs();
    const Size  n_o_d = model.parameters->n_off_diag;

    model.set_dshift_max();//err, probably belongs somewhere else, but we need to compute the max shift for each point

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

        source_c_    (i).resize (width);
        source_n_    (i).resize (width);

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
            //Approach which just accumulates the intensity contributions as the ray is traced

            solve_shortchar_order_0 (model, o, rr);
            solve_shortchar_order_0 (model, o, ar);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            }

            //Approach which first traces a ray

            // const Real dshift_max = get_dshift_max (model, o);
            //
            // nr_   ()[centre] = o;
            // shift_()[centre] = 1.0;
            //
            // first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            // last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            // n_tot_() = (last_()+1) - first_();
            //
            // if (n_tot_() > 1)
            // {
            //     solve_shortchar_order_0_ray_forward (model, o, rr);
            //     solve_shortchar_order_0_ray_backward (model, o, ar);
            //
            //     for (Size f = 0; f < model.parameters->nfreqs(); f++)
            //     {
            //         model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
            //         model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            //     }
            // }
            // else
            // {
            //     for (Size f = 0; f < model.parameters->nfreqs(); f++)
            //     {
            //         model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            //         model.radiation.v(rr,o,f)  = 0.0;
            //     }
            // }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}

/// Solves the radiative transfer equation on the ray using the feautrier solver
// TEMPLATE OPTIONS: (see limitations in assert statements)
//  -IS_SPARSE: whether to use the sparse variant
//  -COMPUTE_UV: whether to compute u, v
//  -COMPUTE_ANIS: whether to compute anisotropic intensity
//  -COMPUTE_LAMBDA: whether to compute the lambda operator
//  -COMPUTE_COOLING_J: whether to compute the mean intensity to compute the cooling rates
template<ApproximationType approx, bool IS_SPARSE, bool COMPUTE_UV, bool COMPUTE_ANIS, bool COMPUTE_LAMBDA, bool COMPUTE_COOLING>
inline void Solver :: solve_feautrier_order_2 (Model& model)
{
    // A few things do currently not make sense, as not everything is implemented in case of sparse/non-sparse solvers
    static_assert(!((!IS_SPARSE)&&(COMPUTE_ANIS)), "Anisotropy not implemented in non-sparse solver.");
    static_assert(!((IS_SPARSE)&&(COMPUTE_UV)), "Computing v not implemented in sparse solver.");
    static_assert(!((COMPUTE_UV)&&(COMPUTE_LAMBDA)), "Computing lambda elements not implemented in uv solver.");
    static_assert(!((COMPUTE_UV)&&(COMPUTE_COOLING)), "Computing cooling not implemented in uv solver");//note: this is due to us not actually computing an angle integrated mean intensity

    ///Intialization for the feautrier solvers
    //Lambda elements need to be cleared before every computation (otherwise we might accidentally sum then with the previous iteration lambda elements)
    if constexpr(COMPUTE_LAMBDA)
    {
        for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}
    }else{}

    if constexpr(COMPUTE_COOLING)
    {
        //reset cooling rates
        for (Size p = 0; p < model.parameters->npoints(); p++)
        {
            model.cooling.cooling_rate[p]=0.0;
        }
    }else{}

    if constexpr(IS_SPARSE)
    {
        //Initialise values for the anisotropic solver
        if constexpr(COMPUTE_ANIS)
        {
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
        }
        else
        {
            // Initialise variables for simple sparse solver
            for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
            {
                lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

                threaded_for (o, model.parameters->npoints(),
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        lspec.J(o,k) = 0.0;
                    }
                })
            }
        }
    }
    else
    {
        //for the non-sparse solvers
        if constexpr(COMPUTE_UV)
        {
            //initialization for the uv solver
            // Allocate memory if not pre-allocated
            if (!model.parameters->store_intensities)
            {
                model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
                model.radiation.v.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
            }
        }
        else
        {
            //initialization for the simplest feautrier solver
            // Allocate memory if not pre-allocated
            if (!model.parameters->store_intensities)
            {
                model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
                model.radiation.J.resize (                            model.parameters->npoints(), model.parameters->nfreqs());
            }

            // Initialise mean intensity
            model.radiation.initialize_J();
        }
    }

    /// Solving the transfer equation for each ray
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];//index of antipod ray

        //c++ being annoying again: constexpr if does introduce local scopes, so I cannot just conditionally declares these variables...
        //But at the same time, it takes almost no time to compute them (and in the base case, the compiler optimizes them away)
        // IDEALLY, something like this would be valid c++ for predefining these variables
        // if constexpr(IS_SPARSE)
        // {
        //     //angular weights for the sparse solvers
        //     const Real wt = model.geometry.rays.weight[rr];//TODO: check if this is the correct weight for the angle
        //     if constexpr(COMPUTE_ANIS)//anisotropic solver needs some extra weights
        //     {
        //         const Vector3D nn = model.geometry.rays.direction[rr];
        //         const Real wt_0    =     inv_sqrt2 * (three * nn.z() * nn.z() - one);
        //         const Real wt_1_Re =        -sqrt3 *          nn.x() * nn.z();
        //         const Real wt_1_Im =        -sqrt3 *          nn.y() * nn.z();
        //         const Real wt_2_Re =  half * sqrt3 * (nn.x() * nn.x() - nn.y() * nn.y());
        //         const Real wt_2_Im =         sqrt3 *          nn.x() * nn.y();
        //     }
        //     else{}
        // }else{}

        //Angular weights for the sparse solvers
        const Real wt      = model.geometry.rays.weight[rr];//TODO: check if this is the correct weight for the angle
        const Vector3D nn  = model.geometry.rays.direction[rr];
        const Real wt_0    =     inv_sqrt2 * (three * nn.z() * nn.z() - one);
        const Real wt_1_Re =        -sqrt3 *          nn.x() * nn.z();
        const Real wt_1_Im =        -sqrt3 *          nn.y() * nn.z();
        const Real wt_2_Re =  half * sqrt3 * (nn.x()* nn.x() - nn.y() * nn.y());
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

                //for the sparse solvers, we iterate over the lines and their quadrature elements
                if constexpr(IS_SPARSE)
                {
                    if (n_tot_() > 1)
                    {
                        for (Size k = 0; k < lspec.linedata.nrad; k++)
                        {
                            // Integrate over the line
                            for (Size z = 0; z < model.parameters->nquads(); z++)
                            {
                                //default feautrier solver is used for sparse solvers
                                solve_feautrier_order_2 <approx> (model, o, lspec.nr_line[o][k][z]);
                                const Real du = lspec.quadrature.weights[z] * wt * Su_()[centre];
                                //For the sparse solvers, we need to compute J (lmean line intensity)
                                lspec.J (o,k) += two * du;
                                if constexpr(COMPUTE_ANIS)
                                {
                                    lspec.J2_0    (o,k) += wt_0    * du;
                                    lspec.J2_1_Re (o,k) += wt_1_Re * du;
                                    lspec.J2_1_Im (o,k) += wt_1_Im * du;
                                    lspec.J2_2_Re (o,k) += wt_2_Re * du;
                                    lspec.J2_2_Im (o,k) += wt_2_Im * du;
                                }
                                else{}

                                //after the solver has done its computations on a single ray, we can compute the lambda values
                                if constexpr(COMPUTE_LAMBDA)
                                {
                                    update_Lambda<approx> (model, rr, lspec.nr_line[o][k][z]);
                                }else{}
                            }
                        }
                    }
                    else
                    {
                        //In case we need a boundary condition for a ray with a single point
                        for (Size k = 0; k < lspec.linedata.nrad; k++)
                        {
                            // Integrate over the line
                            for (Size z = 0; z < model.parameters->nquads(); z++)
                            {
                                const Real du = lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                                //For the sparse solvers, we need to compute J (lmean line intensity)
                                lspec.J (o,k) += two * du;
                                if constexpr(COMPUTE_ANIS)
                                {
                                    lspec.J2_0    (o,k) += wt_0    * du;
                                    lspec.J2_1_Re (o,k) += wt_1_Re * du;
                                    lspec.J2_1_Im (o,k) += wt_1_Im * du;
                                    lspec.J2_2_Re (o,k) += wt_2_Re * du;
                                    lspec.J2_2_Im (o,k) += wt_2_Im * du;
                                }
                                else{}
                            }
                        }
                    }
                }
                //for the non-sparse solvers, we iterate over all frequencies instead
                else
                {
                    if (n_tot_() > 1)
                    {
                        for (Size f = 0; f < model.parameters->nfreqs(); f++)
                        {
                            if constexpr(COMPUTE_UV)
                            {
                                solve_feautrier_order_2_uv <approx> (model, o, f);

                                model.radiation.u(rr,o,f)  = Su_()[centre];
                                model.radiation.v(rr,o,f)  = Sv_()[centre];
                            }
                            else
                            {
                                solve_feautrier_order_2 <approx> (model, o, f);

                                model.radiation.u(rr,o,f)  = Su_()[centre];
                                model.radiation.J(   o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];
                            }

                            //after the solver has done its computations on a single ray, we can compute the lambda values
                            if constexpr(COMPUTE_LAMBDA)
                            {
                                update_Lambda<approx> (model, rr, f);
                            }else{};
                        }
                    }
                    else
                    {
                        for (Size f = 0; f < model.parameters->nfreqs(); f++)
                        {
                            model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                            if constexpr(COMPUTE_UV)
                            {
                                model.radiation.v(rr,o,f)  = 0.0;
                            }
                            else
                            {
                                model.radiation.J(o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
                            }
                        }
                    }
                }
            })
        }
    }

    /// Post-processing

    if constexpr(COMPUTE_COOLING)
    {
        // J contains right the correct amount of information to compute line cooling
        if constexpr(IS_SPARSE)
        {
            threaded_for (p, model.parameters->npoints(),
            {
                for (Size l = 0; l < model.parameters->nlspecs(); l++)
                {
                    const LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        const Size lid = model.lines.line_index(l, k);
                        //line cooling rate maybe gives some useful diagonstics, but can be implemented later on.
                        //line_cooling_rate=same formula
                        model.cooling.cooling_rate[p]+=model.lines.emissivity(p, lid)-lspec.J(p,k)*model.lines.opacity(p, lid);
                    }
                }
            });
        }
        else
        {
            threaded_for (p, model.parameters->npoints(),
            {
                for (Size l = 0; l < model.parameters->nlspecs(); l++)
                {
                    const LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        //quicky compute mean line intensity by integrating over the quadrature
                        //this is much cheaper than evaluating the profile functions all the time
                        Real Jlin=0.0;
                        const Size lid = model.lines.line_index(l, k);
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            Jlin += lspec.quadrature.weights[z] * model.radiation.J(p, lspec.nr_line[p][k][z]);
                        }
                        //line cooling rate maybe gives useful diagonstic, but can be implemented later on.
                        //line_cooling_rate=same formula
                        model.cooling.cooling_rate[p]+=model.lines.emissivity(p, lid)-Jlin*model.lines.opacity(p, lid);
                    }
                }
            });
        }
        //cooling per line is just some manner of diagnostic; total cooling rate might be more interesting

    }else{};

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

    model.   S_ray.resize (n_tot_(), model.parameters->nfreqs());
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

//Because of the new method for computing the optical depth, adding extra frequency points for counteracting the large doppler shift is no longer necessary
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

//Because of the new method for computing the optical depth, adding extra frequency points for counteracting the large doppler shift is no longer necessary
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

    nr   [id1] = nxt;
    shift[id1] = shift_nxt;
    dZ   [id2] = dZ_loc;

    id1 += increment;
    id2 += increment;
}

///  Gaussian line profile function
///    @param[in] width : profile width
///    @param[in] diff  : frequency difference with line centre
///    @return profile function evaluated with this frequency difference
////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: gaussian (const Real inverse_width, const Real diff) const
{
    const Real sqrt_exp = inverse_width * diff;

    return inverse_width * INVERSE_SQRT_PI * expf (-sqrt_exp*sqrt_exp);
}


///  Planck function
///    @param[in] temp : temperature of the corresponding black body
///    @param[in] freq : frequency at which to evaluate the function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: planck (const Real temp, const Real freq) const
{
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1f (HH_OVER_KB*freq/temp);
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
    chi = model.parameters->min_opacity;

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
    chi = prof * model.lines.opacity   (p, l) + model.parameters->min_opacity;
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
    const Size   r)
{
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& source_c= source_c_();
    Vector<Real>& source_n= source_n_();

    Vector<Real>& tau = tau_();


    double  Z = 0.0;   // distance along ray
    double dZ = 0.0;   // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);
    Real term_c, term_n, dtau;
    bool compute_curr_opacity, prev_compute_curr_opacity;

    if (model.geometry.valid_point (nxt))
    {
        double shift_c = 1.0;
        double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);
            const Size l    = model.radiation.frequencies.corresponding_line[f];

            compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

            compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
            dtau = std::max(model.parameters->min_dtau, dtau);

            //proper implementation of 2nd order shortchar (not yet times reducing factor of exp(-tau))
            // model.radiation.I(r,o,f) = term_c * (expm1(-dtau)+dtau) / dtau
            //                          + term_n * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau;
            //Rewrite, trying to use less exponentials
            const Real factor = expm1f(-dtau)/dtau;

            model.radiation.I(r,o,f) = factor*(source_c[f]-source_n[f]*(1.0+dtau))
                                     + source_c[f] - source_n[f];
            tau[f] = dtau;

            //Compute local lambda operator // TODO: check this implementatino when actually using this as a solver for NLTE intensities
            // const Size k = model.radiation.frequencies.corresponding_k_for_tran[f];   // index of transition
            // const Size z = model.radiation.frequencies.corresponding_z_for_line[f];   // index of quadrature point
            // const Real w_ang = model.geometry.rays.weight[r];
            //
            // LineProducingSpecies &lspec = model.lines.lineProducingSpecies[l];
            //
            // const Real freq_line = lspec.linedata.frequency[k];
            // const Real invr_mass = lspec.linedata.inverse_mass;
            // const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;
            //
            // Real inverse_chi=1.0/chi_c[f];
            // Real phi = model.thermodynamics.profile(invr_mass, o, freq_line, freq);
            // Real L   = constante * freq * phi * (factor + 1.0) * inverse_chi;
            // lspec.lambda.add_element(o, k, o, L);
        }

        //For all frequencies, we need to use the same method for computing the optical depth
        // bool prev_compute_curr_opacity=compute_curr_opacity;//technically, we could also keep this bool individually for every frequency
        prev_compute_curr_opacity=compute_curr_opacity;//technically, we could also keep this bool individually for every frequency

        while (model.geometry.not_on_boundary (nxt))
        {
            crt     = nxt;
            shift_c = shift_n;

            model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_n);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                source_c[f]=source_n[f];
                chi_c[f]=chi_n[f];
                const Real freq = model.radiation.frequencies.nu(o, f);
                const Size l    = model.radiation.frequencies.corresponding_line[f];

                compute_curr_opacity=prev_compute_curr_opacity;

                compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
                dtau = std::max(model.parameters->min_dtau, dtau);

                //proper implementation of 2nd order shortchar (not yet times reducing factor of exp(-tau))
                // model.radiation.I(r,o,f) += expf(-tau[f]) *
                //                          ( term_c * (expm1(-dtau)+dtau) / dtau
                //                          + term_n * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau);
                //Rewrite, trying to use less exponentials
                model.radiation.I(r,o,f) += expf(-tau[f]) *
                                           (expm1f(-dtau)/dtau*(source_c[f]-source_n[f]*(1.0+dtau))
                                           + source_c[f] - source_n[f]);
                //TODO: check order of addition, as we might be starting with the largest contributions, before adding the smaller ones...
                tau[f] += dtau;
            }

            //save setting for use for all frequencies for the next interval
            prev_compute_curr_opacity=compute_curr_opacity;
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


// accel inline void Solver :: solve_shortchar_order_0_ray_forward (
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
//     // Size nxt = model.geometry.get_next (o, r, o, Z, dZ);
//     Real term_c, term_n, dtau;
//     bool compute_curr_opacity, prev_compute_curr_opacity;
//     prev_compute_curr_opacity=true;
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq = model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model, nr[first_()], freq*shift[first_()]);
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
//         for (Size f = 0; f < model.parameters->nfreqs(); f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq = model.radiation.frequencies.nu(o, f);
//             const Size l    = model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity = prev_compute_curr_opacity; // for the first point, we need to compute both the curr and next opacity (and source)
//             // compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ[idx-1], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
//             dtau = std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n * (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] * (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] * (expm1(-dtau)+dtau) )/ dtau;
//         }
//         //save setting for use for all frequencies for the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
//     }
// }
//
//
// accel inline void Solver :: solve_shortchar_order_0_ray_backward (
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
//     prev_compute_curr_opacity=true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq = model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model, nr[last_()], freq*(shift[last_()]));
//     }
//     double shift_c, shift_n;
//     //interate until we reach the middle point
//     for (Size indexp1=last_(); indexp1>=centre+1; indexp1--)
//     {
//         crt = nr[indexp1];
//         nxt = nr[indexp1-1];
//         shift_c = shift[indexp1];
//         shift_n = shift[indexp1-1];
//
//         for (Size f = 0; f < model.parameters->nfreqs(); f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq = model.radiation.frequencies.nu(o, f);
//             const Size l    = model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity = prev_compute_curr_opacity; // for the first point, we need to compute both the curr and next opacity (and source)
//             // compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ[indexp1-1], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
//             dtau = std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n * (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] * (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] * (expm1(-dtau)+dtau) )/ dtau;
//         }
//
//         //save setting for use for all frequencies for the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
//     }
// }


template<ApproximationType approx>
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
        // Vector<Real  >& inverse_chi = inverse_chi_();

        const Real w_ang = two * model.geometry.rays.weight[rr];

        const Size l = freqs.corresponding_l_for_spec[f];   // index of species
        const Size k = freqs.corresponding_k_for_tran[f];   // index of transition
        const Size z = freqs.corresponding_z_for_line[f];   // index of quadrature point

        LineProducingSpecies &lspec = model.lines.lineProducingSpecies[l];

        const Real freq_line = lspec.linedata.frequency[k];
        const Real invr_mass = lspec.linedata.inverse_mass;
        const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;
        //TODO: approximation if we do not have overlapping lines
        // Real inverse_opacity = 1.0/model.lines.opacity (nr[centre], k); //Inverse line opacity; includes 1/HH_OVER_FOUR_PI
        // Real L   = constante * L_diag[centre] * inverse_opacity;
        Real eta, chi;//eta is dummy var
        get_eta_and_chi <approx>(model, nr[centre], k, freq_line, eta, chi);
        Real inverse_chi=1.0/chi;

        Real frq = freqs.nu(nr[centre], f) * shift[centre];
        Real phi = thermodyn.profile(invr_mass, nr[centre], freq_line, frq);
        Real L   = constante * frq * phi * L_diag[centre] * inverse_chi;

        lspec.lambda.add_element(nr[centre], k, nr[centre], L);

        for (long m = 0; (m < n_off_diag) && (m+1 < n_tot); m++)
        {
            if (centre >= first+m+1) // centre-m-1 >= first
            {
                const long n = centre-m-1;

                //TODO: approximation if we do not have overlapping lines // also check which opacity we need (centre?)
                // inverse_opacity = 1.0/model.lines.opacity (nr[n], k); //Inverse line opacity; no longer 1/HH_OVER_FOUR_PI
                // L   = constante * L_lower(m,n) * inverse_opacity;
                get_eta_and_chi <approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi=1.0/chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_lower(m,n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }

            if (centre+m+1 <= last) // centre+m+1 < last
            {
                const long n = centre+m+1;

                //TODO: approximation if we do not have overlapping lines
                // inverse_opacity = 1.0/model.lines.opacity (nr[n], k); //Inverse line opacity; no includes 1/HH_OVER_FOUR_PI
                // L   = constante * L_upper(m,n) * inverse_opacity;
                get_eta_and_chi <approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi=1.0/chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_upper(m,n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }
        }
    }
}

///  Computer for the optical depth and source function when computing using the formal line integration
///  In case of low velocity differences, almost two times slower
///  For high velocity increments however, this does not need any extra interpolation points (-> way faster) (and some extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its limit value)
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
template<>
inline void Solver :: compute_S_dtau_line_integrated <OneLine> (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext)
{
    dtau=compute_dtau_single_line(model, currpoint, nextpoint, lineidx, currfreq, nextfreq, dZ);
    Scurr=model.lines.emissivity(currpoint, lineidx)/model.lines.opacity(currpoint, lineidx);//current source
    Snext=model.lines.emissivity(nextpoint, lineidx)/model.lines.opacity(nextpoint, lineidx);//next source
    //note: due to interaction with dtau when computing all sources individually, we do need to recompute Scurr and Snext for all position increments
}

///  Computer for the optical depth and source function when computing using the formal line integration
///  In case of low velocity differences, almost two times slower
///  For high velocity increments however, this does not need any extra interpolation points (-> way faster) (and some extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its limit value)
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
template<>
inline void Solver :: compute_S_dtau_line_integrated <None> (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext)
{
    Real sum_dtau=0.0;
    Real sum_dtau_times_Scurr=0.0;
    Real sum_dtau_times_Snext=0.0;
    for (Size l=0; l<model.parameters->nlines(); l++)
    {
        Real line_dtau=compute_dtau_single_line(model, currpoint, nextpoint, l, currfreq, nextfreq, dZ);
        Real line_Scurr=model.lines.emissivity(currpoint, l)/model.lines.opacity(currpoint, l);//current source
        Real line_Snext=model.lines.emissivity(nextpoint, l)/model.lines.opacity(nextpoint, l);//next source
        sum_dtau+=line_dtau;
        sum_dtau_times_Scurr+=line_dtau*line_Scurr;
        sum_dtau_times_Snext+=line_dtau*line_Snext;
    }
    dtau=sum_dtau;
    Scurr=sum_dtau_times_Scurr/sum_dtau;
    Snext=sum_dtau_times_Snext/sum_dtau;
    //note: due to interaction with dtau when computing all sources individually, we do need to recompute Scurr and Snext for all position increments
}

/// Computes the opacity and optical depth in a hybrid manner
///    @param[in/out] compute_curr_opacity: for deciding whether we need to compute the current opacity when using the trapezoidal rule
///    @param[in] currpoint : index of current point
///    @param[in] nextpoint : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] currshift : shift at curr point
///    @param[in] nextshift : shift at next point
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
template<ApproximationType approx>
accel inline void Solver :: compute_source_dtau (Model& model, Size currpoint, Size nextpoint, Size line, Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ, bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr, Real& Snext)
{
    //deciding which optical depth computation to use, depending on the doppler shift
    const double dshift     = next_shift-curr_shift;//shift[first+1]-shift[first];TODO
    const double dshift_abs = fabs (dshift);
    const double dshift_max = std::min(model.dshift_max[currpoint],model.dshift_max[nextpoint]);
    const bool using_large_shift = (dshift_abs > dshift_max);

    //fancy computation for large doppler shifts
    if (using_large_shift)
    {
        compute_curr_opacity=true;
        compute_S_dtau_line_integrated <approx> (model, currpoint, nextpoint, line, curr_freq, next_freq, dZ, dtaunext, Scurr, Snext);
    }
    else
    {
        //default computation using trapezoidal rule
        if (compute_curr_opacity)//fancy computation does not compute the current opacity, so we might need to recompute it here
        {
            compute_curr_opacity=false;
            Real eta_c=0.0;//current emissivity
            //also get previous opacity (emissivity does not matter)
            get_eta_and_chi <approx> (model, currpoint, line, curr_freq, eta_c, chicurr);
            Scurr=eta_c/chicurr;//might as well compute the source function too
        }

        Real eta_n=0.0;
        // Get new radiative properties
        get_eta_and_chi <approx> (model, nextpoint, line, next_freq, eta_n, chinext);

        Snext = eta_n / chinext;

        dtaunext = half * (chicurr + chinext) * dZ;
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

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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


///   Computes the optical depth assuming only a single line exists.
///    @param[in] curridx: index of the current point
///    @param[in] nextidx: index of the next point
///    @param[in] lineidx: index of the line for which to compute the optical depth
///    @param[in] curr_freq: current frequency (in comoving frame)
///    @param[in] next_freq: next frequency (in comoving frame)
///    @param[in] dz: distance increment
inline Real Solver :: compute_dtau_single_line(Model& model, Size curridx, Size nextidx, Size lineidx, Real curr_freq, Real next_freq, Real dz)
{
    const Real linefreq=model.lines.line[lineidx];
    const Real average_inverse_line_width=(model.lines.inverse_width(curridx, lineidx)+model.lines.inverse_width(nextidx, lineidx))/2.0;

    //opacity is stored divided by the linefreq, so multiply by it
    const Real curr_line_opacity=linefreq*model.lines.opacity(curridx, lineidx);
    const Real next_line_opacity=linefreq*model.lines.opacity(nextidx, lineidx);

    //if frequencies are equal, division by zero (due to the optical depth formula) happens if we were not to use this branch
    if (curr_freq==next_freq)
    {
        //doing the default computation instead (no shifting)
        const Real diff = curr_freq - model.lines.line[lineidx];//curr_freq==next_freq, so choice is arbitrary
        const Real prof = gaussian(average_inverse_line_width, diff);
        const Real average_opacity = (curr_line_opacity+next_line_opacity)/2.0;

        return dz * (prof * average_opacity + model.parameters->min_opacity);
    }

    //We assume a linear interpolation of these dimensionless frequency positions
    //We will also assume the line width to be somewhat constant, replacing the values with the averages
    const Real next_pos=(linefreq-next_freq)*average_inverse_line_width;
    const Real curr_pos=(linefreq-curr_freq)*average_inverse_line_width;

    //In this way, the diff_pos can be computed quite simple, and we do not have a discrepancy between the interpolation and the bounds
    const Real diff_pos=next_pos-curr_pos;

    /// the more correct approach, taking into account also the line opacity change; however, it does not make too much of a difference in the actual result and is quite a bit slower
    // const Real delta_opacity=(next_line_opacity-curr_line_opacity);
    // const Real deltanu=-next_freq+curr_freq;//differences in curr freqs; +-1 due to shift being defined in the other direction
    //
    // //note: opacity can also be extrapolated; however the correction term (expterm) accounts for that
    // const Real interp_opacity=curr_line_opacity+delta_opacity*(curr_freq-linefreq)/deltanu;
    //
    // //This term is a constant term, giving the usual ... as if the opacity were static
    // const Real erfterm=interp_opacity/diff_pos/2.0*(std::erf(next_pos)-std::erf(curr_pos));
    // //This term corrects for the fact that the opacity between points changes
    // const Real expterm=delta_opacity/2.0*INVERSE_SQRT_PI/diff_pos/diff_pos*(std::exp(-curr_pos*curr_pos)-std::exp(-next_pos*next_pos));
    // return dz*std::max(average_inverse_line_width*(erfterm+expterm), model.parameters->min_opacity);

    //If we instead use an average opacity, the computation is quite a bit faster
    const Real average_opacity=(next_line_opacity+curr_line_opacity)/2.0;
    const Real erfterm=average_opacity/diff_pos/2.0*(std::erff(next_pos)-std::erff(curr_pos));
    //correcting to bound opacity from below to the minimum opacity (assumes positive opacities occuring in the model)
    return dz*(average_inverse_line_width*erfterm+model.parameters->min_opacity);

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

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<None>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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

        compute_source_dtau<None>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<None>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    //err, source function might be slightly different when looking at it from curr and next point
    // this is due to the weighting by the line optical depths; this might not be saved TODO think whether this is correct
    model.   S_ray(0, f) =  term_c;
    model.   S_ray(1, f) =  term_n;

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

        compute_source_dtau<None>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        model. S_ray(n+1-first, f) =  term_n;
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

    //a bit fewer vectors are defined than in the default case

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

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
