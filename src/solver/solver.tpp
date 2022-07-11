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

                            update_Lambda (model, rr, lspec.nr_line[o][k][z]);
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
    const Real      dshift_max,
    const int       increment,
          Size      id1,
          Size      id2 )
{
    Real  Z = 0.0;   // distance from origin (o)
    Real dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        Size       crt = o;
        Real shift_crt = geometry.get_shift <frame> (o, r, crt, 0.0);
        Real shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z  );

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


accel inline void Solver :: set_data (
    const Size  crt,
    const Size  nxt,
    const Real  shift_crt,
    const Real  shift_nxt,
    const Real  dZ_loc,
    const Real  dshift_max,
    const int   increment,
          Size& id1,
          Size& id2 )
{
    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

    const Real dshift     = shift_nxt - shift_crt;
    const Real dshift_abs = fabs (dshift);

    if (dshift_abs > dshift_max) // If velocity gradient is not well-sampled enough
    {
        // Interpolate velocity gradient field
        const Size      n_interpl = dshift_abs / dshift_max + 1;
        const Size half_n_interpl = 0.5 * n_interpl;
        const Real     dZ_interpl =     dZ_loc / n_interpl;
        const Real dshift_interpl =     dshift / n_interpl;

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
accel inline Real trap (const Real x_crt, const Real x_nxt, const Real dZ)
{
    return half * (x_crt + x_nxt) * dZ;
}





accel inline void Solver :: solve_shortchar_order_0 (
          Model& model,
    const Size o,
    const Size r,
    const Real dshift_max)
{
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& tau = tau_();


    Real  Z = 0.0;   // distance along ray
    Real dZ = 0.0;   // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point (nxt))
    {
        Real shift_c = 1.0;
        Real shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);

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

        Vector<Size>& nr          = nr_         ();
        Vector<Real>& shift       = shift_      ();
        Vector<Real>& L_diag      = L_diag_     ();
        Matrix<Real>& L_upper     = L_upper_    ();
        Matrix<Real>& L_lower     = L_lower_    ();
        Vector<Real>& inverse_chi = inverse_chi_();

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

    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

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
    // const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);
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

    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

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

    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

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

    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

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

    Vector<Real>& dZ    = dZ_   ();
    Vector<Size>& nr    = nr_   ();
    Vector<Real>& shift = shift_();

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


inline void Solver :: set_eta_and_chi (Model& model, const Size rr) const
{
    model.eta.resize (model.parameters->npoints(), model.parameters->nfreqs());
    model.chi.resize (model.parameters->npoints(), model.parameters->nfreqs());

    for (Size p = 0; p < model.parameters->npoints(); p++)
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            // Extract the Doppler shift
            const Real shift = model.geometry.get_shift <Rest> (0, rr, p, 0.0);
            const Real freq  = model.radiation.frequencies.nu(0, f);
            const Size l     = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi <None> (model, p, l, freq*shift, model.eta(p,f), model.chi(p,f));
        }
    }
}


inline void Solver :: set_boundary_condition (Model& model) const
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

    Real  Z = 0.0;   // distance from origin (o)
    Real dZ = 0.0;   // last increment in Z

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


///  PPORTAL imager
///////////////////
inline void Solver :: PORTAL_image (Model& model, const Size rr, const Size l)
{
    Image image = Image(model.geometry, PolarizedIntensity, rr);

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters->npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            const Size first = first_();
            const Size last  = last_ ();
            const Size n_tot = n_tot_();

            Vector<Real>& dZ    = dZ_   ();
            Vector<Size>& nr    = nr_   ();
            Vector<Real>& shift = shift_();

            image.I_p(o,f) = 0.5*boundary_intensity (model, nr[first], freq*shift[first]);
            image.I_o(o,f) = 0.5*boundary_intensity (model, nr[first], freq*shift[first]);
            image.U  (o,f) = 0.0;

            for (Size n = first; n <= last; n++)
            {
                const Real k_abs_0 = model.k_abs_0[nr[n]];
                const Real k_stm_0 = model.k_stm_0[nr[n]];
                const Real k_abs_2 = model.k_abs_2[nr[n]];
                const Real k_stm_2 = model.k_stm_2[nr[n]];

                const Real cos_t = model.geometry.rays.direction[rr].dot(model.b[nr[n]]);
                const Real fac   = inv_sqrt2 * (three * cos_t * cos_t - two);

                const Real k_abs_p = k_abs_0 + k_abs_2 * fac;
                const Real k_stm_p = k_stm_0 + k_stm_2 * fac;

                const Real k_abs_o = k_abs_0 + k_abs_2 * inv_sqrt2;
                const Real k_stm_o = k_stm_0 + k_stm_2 * inv_sqrt2;

                const Real k_p = k_abs_p - k_stm_p;
                const Real k_o = k_abs_o - k_stm_o;

                const Real freq_ij = model.lines.line[l];
                const Real pre     = HH_OVER_CC_SQUARED * freq_ij * freq_ij * freq_ij;

                const Real S_p = pre * k_stm_p / k_p;
                const Real S_o = pre * k_stm_o / k_o;

                const Real diff = freq*shift[n] - freq_ij;
                const Real prof = gaussian (model.lines.inverse_width(nr[n], l), diff);

                const Real ds_nu = prof * dZ[n];

                const Real dtau_p = k_p * ds_nu;
                const Real dtau_o = k_o * ds_nu;

                image.I_p(o,f) += (image.I_p(o,f) - S_p) * expm1(-dtau_p);
                image.I_o(o,f) += (image.I_o(o,f) - S_o) * expm1(-dtau_o);
                image.U  (o,f) *= exp(-half*(dtau_p + dtau_o));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back (image);
}
