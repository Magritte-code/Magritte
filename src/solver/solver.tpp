///  Getter for the maximum allowed shift value determined by the smallest line
///    @param[in] o : number of point under consideration
///    @retrun maximum allowed shift value determined by the smallest line
///////////////////////////////////////////////////////////////////////////////
// accel inline Real Solver :: get_dshift_max (const Model& model, const Size o) const
// {
//
//     Real dshift_max = std::numeric_limits<Real>::max();
//
//     for (const LineProducingSpecies &lspec : lines.lineProducingSpecies)
//     {
//         const Real inverse_mass   = lspec.linedata.inverse_mass;
//         const Real new_dshift_max = parameters.max_width_fraction
//                                     * thermodynamics.profile_width (inverse_mass, o);
//
//         if (dshift_max > new_dshift_max)
//         {
//             dshift_max = new_dshift_max;
//         }
//     }
//
//     return dshift_max;
// }


inline void Solver :: trace (Model& model)
{
    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "rr = " << rr << endl;

        accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        {
            // const Real dshift_max = get_dshift_max (o);
            const Real dshift_max = 1.0e+99;

            model.geometry.lengths[model.parameters.npoints()*rr+o] =
                trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, +1, centre+1, centre+1) + 1
              - trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, -1, centre+1, centre  );
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec ();
}


inline void Solver :: solve (Model& model)
{
    model.radiation.initialize_J();

    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        // for (Size o = 0; o < model.parameters.npoints(); o++)
        accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        {
            const Real dshift_max = 1.0e+99;

            solve_0th_order_short_charateristics (model, o, rr, dshift_max);
            solve_0th_order_short_charateristics (model, o, ar, dshift_max);

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


inline void Solver :: solve_2nd_order_Feautrier (Model& model)
{
    model.radiation.initialize_J();

    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        // for (Size o = 0; o < model.parameters.npoints(); o++)
        accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        {
            const Real dshift_max = 1.0e+99;

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

           // cout << "   rr = " << rr << "   o = " << o << "   first = " << first_() << "   last = " << last_() << "   n_tot = " << n_tot_() << endl;
           // for (Size i = 0; i < length; i++)
           // {
               // if (i >= first_() && i <= last_())
               // {
                   // cout << " ---- " << i << ":  " << nr_()[i] << "  dZ = " << dZ_()[i] << "   << shift = " << shift_()[i] << endl;
               // }
               // else
               // {
                   // cout << "  " << i << ":  " << nr_()[i] << "  dZ = " << dZ_()[i] << "   << shift = " << shift_()[i] << endl;
               // }
           // }

            //cout << "f = " << first_() << "   l = " << last_() << "   n = " << n_tot_() << endl;

            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                solve_2nd_order_Feautrier (model, o, rr, ar, f);
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
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
    double  Z = 0.0;                  // distance from origin (o)
    double dZ = 0.0;                  // last increment in Z


    // cout << " ---> o = " << o << "  r = " << r << endl;

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    // cout << "          nxt = " << nxt << endl;


    if (geometry.valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt);

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        // cout << "          set data." << endl;

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
                  nxt = geometry.get_next (o, r, nxt, Z, dZ);
            if (nxt == geometry.parameters.npoints()) printf ("ERROR (nxt < 0): o = %ld, crt = %ld, ray = %ld", o, crt, r);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt);

            // cout << "          nxt = " << nxt << endl;

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

            // cout << "          set data." << endl;
        }
    }

    // cout << "          id = " << id << endl;

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
            printf ("ERROR (N_intpl > 10 000) || (dshift_max < 0, probably due to overflow)\n");
        }

        // Assign current cell to first half of interpolation points
        for (Size m = 1; m < half_n_interpl; m++)
        {
            // cout << "id = " << id << endl;

            nr   [id1] = crt;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }

        // Assign next cell to second half of interpolation points
        for (Size m = half_n_interpl; m <= n_interpl; m++)
        {
            // cout << "id = " << id << endl;

            nr   [id1] = nxt;
            shift[id1] = shift_crt + m*dshift_interpl;
            dZ   [id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }
    }

    else
    {
        // cout << "id1 = " << id1 << endl;
        // if (id1 >= length) throw "error";

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

    // cout << "-------------------------------------" << endl;
    // cout << "p = " << p      << endl;
    // cout << "b = " << bdy_id << endl;

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
          Real&  chi )
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
    return half * (x_crt + x_nxt) * dZ;
}





accel inline void Solver :: solve_0th_order_short_charateristics (
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

    // printf("---------- Start: r = %ld,  o = %ld,   nxt = %ld\n", r, o, nxt);

    if (model.geometry.valid_point (nxt))
    {
        double shift_c = 1.0;
        double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt);
        // printf("o = %ld,   nxt = %ld,   shift_nxt = %le\n",o, nxt, shift_nxt);

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            get_eta_and_chi (model, crt, freq,         eta_c[f], chi_c[f]);
            get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);

            const Real drho = trap (eta_c[f], eta_n[f], dZ);
            const Real dtau = trap (chi_c[f], chi_n[f], dZ);

            tau[f]                   = dtau;
            model.radiation.I(r,o,f) = drho * expf(-tau[f]);
            // printf("I(o=%ld, c=%ld, n=%ld) = %le\n", o, crt, nxt, model.radiation.I(r,o,f));
        }

        while (model.geometry.not_on_boundary (nxt))
        {
            crt     = nxt;
            shift_c = shift_n;
              eta_c =   eta_n;
              chi_c =   chi_n;

            model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_n);
            // printf("shift_nxt = %le\n", shift_nxt);

            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                const Real freq = model.radiation.frequencies.nu(o, f);

                get_eta_and_chi (model, nxt, freq*shift_n, eta_n[f], chi_n[f]);

                const Real drho = trap (eta_c[f], eta_n[f], dZ);
                const Real dtau = trap (chi_c[f], chi_n[f], dZ);

                tau[f]                   += dtau;
                model.radiation.I(r,o,f) += drho * expf(-tau[f]);
                // printf("I(o=%ld, c=%ld, n=%ld) = %le\n", o, crt, nxt, model.radiation.I(r,o,f));
            }
        }
        // printf("shift_nxt = %le\n", shift_nxt);

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f) += boundary_intensity(model, nxt, freq*shift_n) * expf(-tau[f]);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
            // printf("I(o=%ld, c=%ld, n=%ld) = %le\n", o, crt, nxt, model.radiation.I(r,o,f));
            // printf("-------- bc(%ld, %le) = %le\n", nxt, freq*shift_nxt, boundary_intensity(model, nxt, freq*shift_nxt));
        }
    }

    else
    {
        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f)  = boundary_intensity(model, crt, freq);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
            // printf("I(o=%ld, c=%ld, n=%ld) = %le\n", o, crt, nxt, model.radiation.I(r,o,f));
            // printf("-------- bc(%ld, %le) = %le\n", crt, freq, boundary_intensity(model, crt, freq));
        }
    }
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
accel inline void Solver :: solve_2nd_order_Feautrier (
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

    term_c = eta_c / chi_c;
    term_n = eta_n / chi_n;
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
        get_eta_and_chi (model, nr[n], freq*shift[n], eta_n, chi_n);

        term_n = eta_n / chi_n;
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

//            printf("- Solver: first %ld, last %ld, n_tot %ld, r %ld\n", first, last, n_tot[rp], rr);
//            printf("L(n=%ld, f=%ld) = %le\n", centre, f, L_diag[In1]);
        }
        else
        {
            L_diag[centre] = (one + FF[centre-1]) / (Bl_min_Al + Bl*FF[centre-1]);
//            printf("- Solver: first %ld, last %ld, n_tot %ld, r %ld\n", first, last, n_tot[rp], rr);
//            printf("L(n=%ld, f=%ld) = %le\n", centre, f, L_diag[In1]);
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
//            printf("L(n=%ld, f=%ld) = %le\n", n, f, L_diag[In]);
        }

        Su    [first] += Su[first+1] * FI[first];
        L_diag[first]  = (one + GG[first+1]) / (Bf_min_Cf + Bf*GG[first+1]);

//        printf("first G[Ifp1] = %le\n", G[Ifp1]);
//        printf("L(n=%ld, f=%ld) = %le\n", first, f, L_diag[If]);


        for (long n = last-1; n >= first; n--) // use long in reverse loops!
        {
            L_upper(0,n+1) = L_diag[n+1] * FI[n  ];
            L_lower(0,n  ) = L_diag[n  ] * GI[n+1];

//            printf("L_u(0, %ld) = %le\n", In, L_upper[M(0,In)]);
//            printf("L_l(0, %ld) = %le\n", In, L_lower[M(0,In)]);
        }

        for (Size m = 1; (m < n_off_diag) && (m < n_tot-1); m++)
        {
            for (long n = last-1-m; n >= first; n--) // use long in reverse loops!
            {
                L_upper(m,n+m+1) = L_upper(m-1,n+m+1) * FI[n    ];
                L_lower(m,n    ) = L_lower(m-1,n    ) * GI[n+m+1];

//                printf("L_u(%ld, %ld)[%ld] = %le\n", m, In, M(m,In), L_upper[M(m,In)]);
//                printf("L_l(%ld, %ld)[%ld] = %le\n", m, In, M(m,In), L_lower[M(m,In)]);
            }
        }

//        for (long n = last-1; n >= first; n--) // use long in reverse loops!
//        {
//            const Size Inp1 = I(n+1, w);
//            const Size In   = I(n,   w);
//
//            printf("check L_u(0, %ld)[%ld] = %le\n", In, M(0,In), L_upper[M(0,In)]);
//            printf("check L_l(0, %ld)[%ld] = %le\n", In, M(0,In), L_lower[M(0,In)]);
//        }
    }



    model.radiation.u(rr,o,f)  = Su[centre];
    model.radiation.J(   o,f) += Su[centre] * two * model.geometry.rays.weight[rr];


}
