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

        accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
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


// inline void Solver :: trace (Model& model)
// {
//     for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
//     {
//         const Size ar = model.geometry.rays.antipod[rr];
//
//         cout << "rr = " << rr << endl;
//
//         accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
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

        //accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        for (Size o = 0; o < model.parameters.npoints(); o++)
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
        }

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}


inline void Solver :: solve_feautrier_order_2 (Model& model)
{

    // Clear the approximate LAMBDA operator
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    model.radiation.initialize_J();

    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        //accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        for (Size o = 0; o < model.parameters.npoints(); o++)
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                //for (Size f = 0; f < model.parameters.nfreqs(); f++)
                accelerated_for (f, model.parameters.nfreqs(), nblocks, nthreads,
                {
                    solve_feautrier_order_2 (model, o, rr, ar, f);

                    model.radiation.u(rr,o,f)  = Su_()[centre];
                    model.radiation.J(   o,f) += Su_()[centre] * TWO * model.geometry.rays.weight[rr];

                    update_Lambda (model, rr, f);
                })
            }
            else
            {
                for (Size f = 0; f < model.parameters.nfreqs(); f++)
                {
                    model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.J(   o,f) += TWO * model.geometry.rays.weight[rr] * model.radiation.u(rr,o,f);
                }
            }
        }

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

    accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = 1.0;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                image_feautrier_order_2 (model, o, rr, ar, f);

                image.I(o,f) = TWO*Su_()[first_()] - boundary_intensity(model, nr_()[first_()], model.radiation.frequencies.nu(o, f));
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


///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap (const Real x_crt, const Real x_nxt, const double dZ)
{
    return HALF * (x_crt + x_nxt) * dZ;
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

        const Real w_ang = TWO * model.geometry.rays.weight[rr];

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
    dtau_n = HALF * (chi_c + chi_n) * dZ[first];

    // Set boundary conditions
    const Real inverse_dtau_f = ONE / dtau_n;

            C[first] = TWO * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = 1.0 / C[first];   // Required for Lambda_diag

    const Real Bf_min_Cf = ONE + TWO * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + TWO * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = HALF * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = ONE / (ONE + FF[first]);


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
        dtau_n = HALF * (chi_c + chi_n) * dZ[n];

        const Real dtau_avg = HALF * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = ONE / inverse_A[n];
        C[n] = ONE / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + ONE) * inverse_C[n];
        FI[n] = ONE / (ONE + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = ONE / dtau_n;

    A[last] = TWO * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = ONE + TWO * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = ONE / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + TWO * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (ONE + FF[last-1]) * denominator;

    if (n_off_diag == 0)
    {
        if (centre < last)
        {
            /// Write economically: G[last] = (B[last] - A[last]) / A[last];
            GG[last] = HALF * Bl_min_Al * dtau_n * dtau_n;
            GP[last] = GG[last] / (ONE + GG[last]);

            for (long n = last-1; n > centre; n--) // use long in reverse loops!
            {
                Su[n] += Su[n+1] * FI[n];

                GG[n] = (C[n] * GP[n+1] + ONE) * inverse_A[n];
                GP[n] = GG[n] / (ONE + GG[n]);
            }

            Su    [centre] += Su[centre+1] * FI[centre];
            L_diag[centre]  = inverse_C[centre] / (FF[centre] + GP[centre+1]);
        }
        else
        {
            L_diag[centre] = (ONE + FF[centre-1]) / (Bl_min_Al + Bl*FF[centre-1]);
        }
    }
    else
    {
        /// Write economically: G[last] = (B[last] - A[last]) / A[last];
        GG[last] = HALF * Bl_min_Al * dtau_n * dtau_n;
        GI[last] = ONE / (ONE + GG[last]);
        GP[last] = GG[last] * GI[last];

        L_diag[last] = (ONE + FF[last-1]) / (Bl_min_Al + Bl*FF[last-1]);

        for (long n = last-1; n > first; n--) // use long in reverse loops!
        {
            Su[n] += Su[n+1] * FI[n];

            GG[n] = (C[n] * GP[n+1] + ONE) * inverse_A[n];
            GI[n] = ONE / (ONE + GG[n]);
            GP[n] = GG[n] * GI[n];

            L_diag[n] = inverse_C[n] / (FF[n] + GP[n+1]);
        }

        Su    [first] += Su[first+1] * FI[first];
        L_diag[first]  = (ONE + GG[first+1]) / (Bf_min_Cf + Bf*GG[first+1]);

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
    dtau_n = HALF * (chi_c + chi_n) * dZ[first];

    // Set boundary conditions
    const Real inverse_dtau_f = ONE / dtau_n;

    C[first] = TWO * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = ONE + TWO * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + TWO * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = HALF * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = ONE / (ONE + FF[first]);

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
        dtau_n = HALF * (chi_c + chi_n) * dZ[n];

        // cout << "dtau = " << dtau_n << endl;

        const Real dtau_avg = HALF * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = ONE / inverse_A[n];
        C[n] = ONE / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        // cout << "inverse_C[" << n << "] = " << inverse_C[n] << endl;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + ONE) * inverse_C[n];
        FI[n] = ONE / (ONE + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = ONE / dtau_n;

    A[last] = TWO * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = ONE + TWO * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = ONE / (Bl * FF[last-1] + Bl_min_Al);

    // cout << "Bl = " << Bl << "   FF[last-1] = " << FF[last-1] << "   Bl_min_Al = " << Bl_min_Al << endl;

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + TWO * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (ONE + FF[last-1]) * denominator;

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
