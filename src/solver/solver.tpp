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
                trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, +1)
              - trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, -1);
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec ();
}


inline void Solver :: solve (Model& model)
{
    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "rr = " << rr << endl;

        accelerated_for (o, model.parameters.npoints(), nblocks, nthreads,
        {
            const Real dshift_max = 1.0e+99;

            solve_0th_order_short_charateristics (model, o, rr, dshift_max);
            solve_0th_order_short_charateristics (model, o, ar, dshift_max);
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec ();
}


template <Frame frame>
accel inline Size Solver :: trace_ray (
    const Geometry& geometry,
    const Size      o,
    const Size      r,
    const Real      dshift_max,
    const int       increment )
{
    Size id = centre;   // index distance from origin
    Real  Z = 0.0;      // distance from origin (o)
    Real dZ = 0.0;      // last increment in Z

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (nxt != geometry.parameters.npoints()) // if we are not going out of mesh
    {
        Size       crt = o;
        Real shift_crt = geometry.get_shift <frame> (o, r, crt);
        Real shift_nxt = geometry.get_shift <frame> (o, r, nxt);

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id);

        while (geometry.boundary.point2boundary[nxt] == geometry.parameters.npoints()) // while nxt not on boundary
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
                  nxt = geometry.get_next (o, r, nxt, Z, dZ);
            if (nxt == geometry.parameters.npoints()) printf ("ERROR (nxt < 0): o = %ld, crt = %ld, ray = %ld", o, crt, r);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt);

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id);
        }
    }

    return id;
}


accel inline void Solver :: set_data (
    const Size  crt,
    const Size  nxt,
    const Real  shift_crt,
    const Real  shift_nxt,
    const Real  dZ_loc,
    const Real  dshift_max,
    const int   increment,
          Size& id )
{
    const Real dshift     = shift_nxt - shift_crt;
    const Real dshift_abs = fabs (dshift);

    if (dshift_abs > dshift_max) // If velocity gradient is not well-sampled enough
    {
        // Interpolate velocity gradient field
        const Size      n_interpl = dshift_abs / dshift_max + Real (1);
        const Size half_n_interpl = Real (0.5) * n_interpl;
        const Real     dZ_interpl =     dZ_loc / n_interpl;
        const Real dshift_interpl =     dshift / n_interpl;

        if (n_interpl > 10000)
        {
            printf ("ERROR (N_intpl > 10 000) || (dshift_max < 0, probably due to overflow)\n");
        }

        // Assign current cell to first half of interpolation points
        for (Size m = 1; m < half_n_interpl; m++)
        {
               nr[id] = crt;
            shift[id] = shift_crt + m*dshift_interpl;
               dZ[id] = dZ_interpl;

            id += increment;
        }

        // Assign next cell to second half of interpolation points
        for (Size m = half_n_interpl; m <= n_interpl; m++)
        {
               nr[id] = nxt;
            shift[id] = shift_crt + m*dshift_interpl;
               dZ[id] = dZ_interpl;

            id += increment;
        }
    }

    else
    {
           nr[id] = nxt;
        shift[id] = shift_nxt;
           dZ[id] = dZ_loc;

           id += increment;
    }
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////

// template <typename Real, typename DataLayout>
// HOST_DEVICE
// inline void Solver<Real, DataLayout> :: solve_2nd_order_Feautrier_non_adaptive (const Size w)
// {
//     // Get indices of the block
//     const Size rp = w / nfreqs_red;   // raypair index
//     const Size f  = w % nfreqs_red;   // frequency index
//
//     // Get frequency value corresponding to the block
//     const Real frequency = frequencies[V(origins[rp], f)];
//
//     // Get indices for first element in the block
//     const Size If   = I(first,   w);
//     const Size Ifp1 = I(first+1, w);
//     const Size Df   = D(rp, first);
//
//     // Get optical properties for first two elements
//     get_eta_and_chi (If,   Df,   frequency);
//     get_eta_and_chi (Ifp1, Df+1, frequency);
//
//     term1[If  ] = eta[If]   / chi[If];
//     term1[Ifp1] = eta[Ifp1] / chi[Ifp1];
//
//     // Get first optical depth increment
//     dtau[If] = 0.5 * (chi[If] + chi[Ifp1]) * dZs[Df];
//
//     // Set boundary conditions
//     const Real inverse_dtau0 = one / dtau[If];
//
//     C[If] = 2.0 * inverse_dtau0 * inverse_dtau0;
//
//     const Real B0_min_C0 = my_fma (2.0, inverse_dtau0, one);
//     const Real B0        = B0_min_C0 + C[If];
//
// //    const Real I_bdy_0 = planck (T_CMB, frequency*shifts[Df]);
//     const Real I_bdy_0 = boundary_intensity (bdy_0[rp], frequency*shifts[Df]);
//
//     Su[If] = term1[If] + 2.0 * I_bdy_0 * inverse_dtau0;
//
//     /// Start elimination step
//     Su[If] = Su[If] / B0;
//
//     /// Write economically: F[first] = (B[first] - C[first]) / C[first];
//     F[If] = 0.5 * B0_min_C0 * dtau[If] * dtau[If];
//     inverse_one_plus_F[If] = one / (one + F[If]);
//
//
//     /// Set body of Feautrier matrix
//     for (Size n = first+1; n < last; n++)
//     {
//         const Size Inm1  = I(n-1, w);
//         const Size In    = I(n,   w);
//         const Size Inp1  = I(n+1, w);
//         const Size Dn    = D(rp,  n);
//
//         // Get new optical properties
//         get_eta_and_chi (Inp1, Dn+1, frequency);
//
//         // Compute term1 at n+1
//         term1[Inp1] = eta[Inp1] / chi[Inp1];
//
//         // Compute the optical depth increment
//         dtau[In] = 0.5 * (chi[In] + chi[Inp1]) * dZs[Dn];
//
//         const Real dtau_avg = 0.5 * (dtau[Inm1] + dtau[In]);
//         inverse_A[In] = dtau_avg * dtau[Inm1];
//         inverse_C[In] = dtau_avg * dtau[In  ];
//
//         A[In] = one / inverse_A[In];
//         C[In] = one / inverse_C[In];
//
//         /// Use the previously stored value of the source function
//         Su[In] = term1[In];
//
//
//         F[In] = my_fma (A[In]*F[Inm1], inverse_one_plus_F[Inm1], one) * inverse_C[In];
//         inverse_one_plus_F[In] = one / (one + F[In]);
//
//         Su[In] = my_fma (A[In], Su[Inm1], Su[In]) * inverse_one_plus_F[In] * inverse_C[In];
//     }
//
//
//     // Get indices for first element in the block
//     const Size Il   = I(last,   w);
//     const Size Ilm1 = I(last-1, w);
//     const Size Dl   = D(rp, last);
//
//     /// Set boundary conditions
//     const Real inverse_dtaud = one / dtau[Ilm1];
//
//     A[Il] = 2.0 * inverse_dtaud * inverse_dtaud;
//
//     const Real Bd_min_Ad = my_fma (2.0, inverse_dtaud, one);
//     const Real Bd        = Bd_min_Ad + A[Il];
//
//     const Real denominator = one / my_fma (Bd, F[Ilm1], Bd_min_Ad);
//
// //    const Real I_bdy_n = planck (T_CMB, frequency*shifts[Dl]);
//     const Real I_bdy_n = boundary_intensity (bdy_n[rp], frequency*shifts[Df]);
//
//     Su[Il] = term1[Il] + 2.0 * I_bdy_n * inverse_dtaud;
//     Su[Il] = my_fma (A[Il], Su[Ilm1], Su[Il]) * (one + F[Ilm1]) * denominator;
//
//
//     if (n_off_diag == 0)
//     {
//         if (n1_min < last)
//         {
//             /// Write economically: G[last] = (B[last] - A[last]) / A[last];
//             G[Il] = 0.5 * Bd_min_Ad * dtau[Ilm1] * dtau[Ilm1];
//             G_over_one_plus_G[Il] = G[Il] / (one + G[Il]);
//
//             for (long n = last-1; n > n1_min; n--) // use long in reverse loops!
//             {
//                 const Size Inp1 = I(n+1, w);
//                 const Size In   = I(n,   w);
//
//                 Su[In] = my_fma(Su[Inp1], inverse_one_plus_F[In], Su[In]);
//
//                 G[In] = my_fma(C[In], G_over_one_plus_G[Inp1], one) * inverse_A[In];
//                 G_over_one_plus_G[In] = G[In] / (one + G[In]);
//             }
//
//             const Size In1   = I(n1_min,   w);
//             const Size In1p1 = I(n1_min+1, w);
//
//             Su[In1] = my_fma (Su[In1p1], inverse_one_plus_F[In1], Su[In1]);
//             L_diag[In1] = inverse_C[In1] / (F[In1] + G_over_one_plus_G[In1p1]);
//
// //            printf("- Solver: first %ld, last %ld, n_tot %ld, r %ld\n", first, last, n_tot[rp], rr);
// //            printf("L(n=%ld, f=%ld) = %le\n", n1_min, f, L_diag[In1]);
//         }
//         else
//         {
//             const Size In1   = I(n1_min,   w);
//             const Size In1m1 = I(n1_min-1, w);
//
//             L_diag[In1] = (one + F[In1m1]) / (Bd_min_Ad + Bd*F[In1m1]);
// //            printf("- Solver: first %ld, last %ld, n_tot %ld, r %ld\n", first, last, n_tot[rp], rr);
// //            printf("L(n=%ld, f=%ld) = %le\n", n1_min, f, L_diag[In1]);
//         }
//     }
//     else
//     {
//         /// Write economically: G[last] = (B[last] - A[last]) / A[last];
//         G[Il] = 0.5 * Bd_min_Ad * dtau[Ilm1] * dtau[Ilm1];
//         inverse_one_plus_G[Il] = one / (one + G[Il]);
//         G_over_one_plus_G[Il] = G[Il] * inverse_one_plus_G[Il];
//
//         L_diag[Il] = (one + F[Ilm1]) / (Bd_min_Ad + Bd*F[Ilm1]);
//
//         for (long n = last-1; n > first; n--) // use long in reverse loops!
//         {
//             const Size Inp1 = I(n+1, w);
//             const Size In   = I(n,   w);
//
//             Su[In] = my_fma(Su[Inp1], inverse_one_plus_F[In], Su[In]);
//
//             G[In] = my_fma(C[In], G_over_one_plus_G[Inp1], one) * inverse_A[In];
//             inverse_one_plus_G[In] = one / (one + G[In]);
//             G_over_one_plus_G[In] = G[In] * inverse_one_plus_G[In];
//
//             L_diag[In] = inverse_C[In] / (F[In] + G_over_one_plus_G[Inp1]);
// //            printf("L(n=%ld, f=%ld) = %le\n", n, f, L_diag[In]);
//         }
//
//         Su[If] = my_fma(Su[Ifp1], inverse_one_plus_F[If], Su[If]);
//         L_diag[If] = (one + G[Ifp1]) / (B0_min_C0 + B0*G[Ifp1]);
//
// //        printf("first G[Ifp1] = %le\n", G[Ifp1]);
// //        printf("L(n=%ld, f=%ld) = %le\n", first, f, L_diag[If]);
//
//
//         for (long n = last-1; n >= first; n--) // use long in reverse loops!
//         {
//             const Size Inp1 = I(n+1, w);
//             const Size In   = I(n,   w);
//
//             L_upper[M(0,Inp1)] = L_diag[Inp1] * inverse_one_plus_F[In  ];
//             L_lower[M(0,In  )] = L_diag[In  ] * inverse_one_plus_G[Inp1];
//
// //            printf("L_u(0, %ld) = %le\n", In, L_upper[M(0,In)]);
// //            printf("L_l(0, %ld) = %le\n", In, L_lower[M(0,In)]);
//         }
//
//         for (Size m = 1; (m < n_off_diag) && (m < n_tot[rp]-1); m++)
//         {
//             for (long n = last-1-m; n >= first; n--) // use long in reverse loops!
//             {
//                 const Size Inp1   = I(n+1,   w);
//                 const Size Inpmp1 = I(n+m+1, w);
//                 const Size Inpmp2 = I(n+m+2, w);
//                 const Size In     = I(n,     w);
//
//                 L_upper[M(m,Inpmp1)] = L_upper[M(m-1,Inpmp1)] * inverse_one_plus_F[In    ];
//                 L_lower[M(m,In    )] = L_lower[M(m-1,In    )] * inverse_one_plus_G[Inpmp1];
//
// //                printf("L_u(%ld, %ld)[%ld] = %le\n", m, In, M(m,In), L_upper[M(m,In)]);
// //                printf("L_l(%ld, %ld)[%ld] = %le\n", m, In, M(m,In), L_lower[M(m,In)]);
//             }
//         }
//
// //        for (long n = last-1; n >= first; n--) // use long in reverse loops!
// //        {
// //            const Size Inp1 = I(n+1, w);
// //            const Size In   = I(n,   w);
// //
// //            printf("check L_u(0, %ld)[%ld] = %le\n", In, M(0,In), L_upper[M(0,In)]);
// //            printf("check L_l(0, %ld)[%ld] = %le\n", In, M(0,In), L_lower[M(0,In)]);
// //        }
//     }
//
// }


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
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1 (HH_OVER_KB*freq/temp);
}


accel inline Real Solver :: boundary_intensity (const Model& model, const Size bdy_id, const Real freq) const
{
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
    chi = 0.0;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters.nlines(); l++)
    {
        const Real diff = freq - model.lines.line[l];
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);
    }
}


///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap (const Real x_crt, const Real x_nxt, const Real dZ)
{
    return (Real) 0.5 * (x_crt + x_nxt) * dZ;
}




accel inline void Solver :: solve_0th_order_short_charateristics (
          Model& model,
    const Size   o,
    const Size   r,
    const Real   dshift_max)
{
    Matrix<Real>& I = model.radiation.I[r];

    Real  Z = 0.0;   // distance along ray
    Real dZ = 0.0;   // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point (nxt))
    {
        Real shift_crt = 1.0;
        Real shift_nxt = model.geometry.get_shift <CoMoving> (o, r, nxt);

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            get_eta_and_chi (model, crt, freq,           eta_crt[f], chi_crt[f]);
            get_eta_and_chi (model, nxt, freq*shift_nxt, eta_nxt[f], chi_nxt[f]);

            drho[f] = trap (eta_crt[f], eta_nxt[f], dZ);
            dtau[f] = trap (chi_crt[f], chi_nxt[f], dZ);
             tau[f] = dtau[f];
             I(o,f) = drho[f] * expf(-tau[f]);
        }

        while (model.geometry.not_on_boundary (nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
              eta_crt =   eta_nxt;
              chi_crt =   chi_nxt;

            model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_nxt);

            for (Size f = 0; f < model.parameters.nfreqs(); f++)
            {
                const Real freq = model.radiation.frequencies.nu(o, f);

                get_eta_and_chi (model, nxt, freq*shift_nxt, eta_nxt[f], chi_nxt[f]);

                drho[f]  = trap (eta_crt[f], eta_nxt[f], dZ);
                dtau[f]  = trap (chi_crt[f], chi_nxt[f], dZ);
                 tau[f] += dtau[f];
                 I(o,f) += drho[f] * expf(-tau[f]);
            }
        }

        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            I(o,f) += boundary_intensity(model, nxt, freq*shift_nxt) * expf(-tau[f]);
        }
    }

    else
    {
        for (Size f = 0; f < model.parameters.nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            I(o,f) = boundary_intensity(model, crt, freq);
        }
    }
}
