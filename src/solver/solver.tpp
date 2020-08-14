void Solver :: trace (const Model& model)
{
    for (size_t rr = 0; rr < model.geometry.rays.get_nrays(); rr++)
    {
        const size_t ar = model.geometry.rays.antipod[rr];

        for (size_t o = 0; o < model.geometry.points.get_npoints(); o++)
        {
            const double dshift_max = 1.0e+99;

            trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, +1);
            trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, -1);
        }
    }
}


template <Frame frame>
inline void Solver :: trace_ray (
    const Geometry& geometry,
    const size_t    o,
    const size_t    r,
    const double    dshift_max,
    const int       increment   )
{
    size_t id = centre;   // index distance from origin
    double  Z = 0.0;      // distance from origin (o)
    double dZ = 0.0;      // last increment in Z

    size_t nxt = geometry.get_next (o, r, o, Z, dZ);

    if (nxt != geometry.get_npoints()) // if we are not going out of mesh
    {
        size_t       crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt);

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id);

        while (!geometry.boundary.is_on_boundary[nxt])
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
                  nxt = geometry.get_next (o, r, nxt, Z, dZ);
            if (nxt == geometry.get_npoints()) printf ("ERROR (nxt < 0): o = %ld, crt = %ld, ray = %ld", o, crt, r);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt);

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id);
        }
    }
}


inline void Solver :: set_data (
    const size_t  crt,
    const size_t  nxt,
    const double  shift_crt,
    const double  shift_nxt,
    const double  dZ_loc,
    const double  dshift_max,
    const int     increment,
          size_t& id           )
{
    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs (dshift);

    if (dshift_abs > dshift_max) // If velocity gradient is not well-sampled enough
    {
        // Interpolate velocity gradient field
        const size_t      n_interpl = dshift_abs / dshift_max + 1;
        const size_t half_n_interpl =        0.5 * n_interpl;
        const double     dZ_interpl =     dZ_loc / n_interpl;
        const double dshift_interpl =     dshift / n_interpl;

        if (n_interpl > 10000)
        {
            printf ("ERROR (N_intpl > 10 000) || (dshift_max < 0, probably due to overflow)\n");
        }

        // Assign current cell to first half of interpolation points
        for (size_t m = 1; m < half_n_interpl; m++)
        {
               nr[id] = crt;
            shift[id] = shift_crt + m*dshift_interpl;
               dZ[id] = dZ_interpl;

            id += increment;
        }

        // Assign next cell to second half of interpolation points
        for (size_t m = half_n_interpl; m <= n_interpl; m++)
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



