inline void Geometry :: set_npoints (const size_t n)
{
    npoints = n;
}


inline size_t Geometry :: get_npoints () const
{
    return npoints;
}


///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      o : number of cell from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
inline size_t Geometry :: get_next (
    const size_t  o,
    const size_t  r,
    const size_t  c,
          double& Z,
          double& dZ               ) const
{
    const size_t     n_nbs = points.    n_neighbors[c];
    const size_t cum_n_nbs = points.cum_n_neighbors[c];

    double dmin = std::numeric_limits<double>::max();   // Initialize to "infinity"
    size_t next = npoints;                              // return npoints when there is no next

    for (size_t i = 0; i < n_nbs; i++)
    {
        const size_t n     = points.neighbors[cum_n_nbs+i];
        const Vector R     = points.position[n] - points.position[o];
        const double Z_new = R.dot(rays.direction[r]);

        if (Z_new > Z)
        {
            const double distance_from_ray2 = R.dot(R) - Z_new*Z_new;

            if (distance_from_ray2 < dmin)
            {
                dmin = distance_from_ray2;
                next = n;
                dZ   = Z_new - Z;   // such that dZ > 0.0
            }
        }
    }

    // Update distance along ray
    Z += dZ;

    return next;
}




///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template <Frame frame>
inline double Geometry :: get_shift (
    const size_t  o,
    const size_t  r,
    const size_t  crt               ) const
{
    // Co-moving frame implementation
    if (frame == CoMoving)
    {
        return 1.0 - (points.velocity[crt] - points.velocity[o]).dot(rays.direction[r]);
    }

    // Rest frame implementation
    if (frame == Rest)
    {
        // In the rest frame the direction of the projected should be fixed
        // We choose to fix it to "up the ray"

        size_t r_correct = r;

        if (r >= nrays/2)
        {
            r_correct = rays.antipod[r];
        }

        return 1.0 - points.velocity[crt].dot(rays.direction[r_correct]);
    }
}


//inline double Geometry :: get_shift <CoMoving> (
//        const size_t  o,
//        const size_t  r,
//        const size_t  crt                      ) const
//{
//        return 1.0 - (points.velocity[crt] - points.velocity[o]).dot(rays.direction[r]);
//}
//
//inline double Geometry :: get_shift <Rest> (
//        const size_t  o,
//        const size_t  r,
//        const size_t  crt                  ) const
//{
//    size_t r_correct = r;
//
//    if (r >= nrays/2)
//    {
//        r_correct = rays.antipod[r];
//    }
//
//    return 1.0 - points.velocity[crt].dot(rays.direction[r_correct]);
//}


inline size_t Geometry :: get_n_interpl (
        const double  shift_crt,
        const double  shift_nxt,
        const double  dshift_max      ) const
{
    const double dshift_abs = fabs (shift_nxt - shift_crt);

    if (dshift_abs > dshift_max) {return dshift_abs/dshift_max + 1;}
    else                         {return 1;                        }
}



template <Frame frame>
inline size_t Geometry :: get_ray_length (
        const size_t    o,
        const size_t    r,
        const double    dshift_max     ) const
{
    size_t  l = 0;        // ray length
    double  Z = 0.0;      // distance from origin (o)
    double dZ = 0.0;      // last increment in Z

    size_t nxt = get_next (o, r, o, Z, dZ);

    if (nxt != get_npoints()) // if we are not going out of mesh
    {
        size_t       crt = o;
        double shift_crt = get_shift <frame> (o, r, crt);
        double shift_nxt = get_shift <frame> (o, r, nxt);

        l += get_n_interpl (shift_crt, shift_nxt, dshift_max);

        while (!boundary.is_on_boundary[nxt])
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
                  nxt = get_next (o, r, nxt, Z, dZ);
            if (nxt == get_npoints()) {
                printf ("ERROR (nxt < 0): o = %ld, crt = %ld, ray = %ld\n", o, crt, r);

                throw "PROBLEM";};
            shift_nxt = get_shift <frame> (o, r, nxt);

            l += get_n_interpl (shift_crt, shift_nxt, dshift_max);
        }
    }

    return l;
}


inline Long2 Geometry :: get_ray_lengths () const
{
    const size_t hnrays  = rays  .get_nrays  ()/2;
    const size_t npoints = points.get_npoints();

    Long2 lengths (hnrays, Long1 (npoints, 0));

    for (size_t rr = 0; rr < hnrays; rr++)
    {
        cout << "rr = " << rr << endl;

        const size_t ar = rays.antipod[rr];

//        for (size_t o = 0; o < npoints; o++)
        threaded_for (o, npoints,
        {
            const double dshift_max = 1.0e+99;

            lengths[rr][o] += get_ray_length <CoMoving> (o, rr, dshift_max);
            lengths[rr][o] += get_ray_length <CoMoving> (o, ar, dshift_max);
        })
    }

    return lengths;
}
