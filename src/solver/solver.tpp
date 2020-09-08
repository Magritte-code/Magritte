inline void Solver :: trace (Model& model)
{
    for (Size rr = 0; rr < model.parameters.hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "rr = " << rr << endl;

        accelerated_for (o, 10000, nblocks, nthreads,
        {
            const Real dshift_max = 1.0e+99;

            model.geometry.lengths[model.parameters.npoints()*rr+o] =
                trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, +1)
              - trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, -1);
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
        const Real      n_interpl = dshift_abs / dshift_max + Real (1);
        const Real half_n_interpl = Real (0.5) * n_interpl;
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
