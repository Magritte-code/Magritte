#include <limits>
#include <cmath>


accel inline void Geometry :: set_npoints (const Size n)
{
    npoints = n;
}


accel inline Size Geometry :: get_npoints () const
{
    return npoints;
}


accel inline void Geometry :: set_nrays (const Size n)
{
    nrays = n;
}


accel inline Size Geometry :: get_nrays () const
{
    return nrays;
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
accel inline Size Geometry :: get_next (
    const Size  o,
    const Size  r,
    const Size  c,
          Real& Z,
          Real& dZ                     ) const
{
    const Size     n_nbs = points.    n_neighbors[c];
    const Size cum_n_nbs = points.cum_n_neighbors[c];

    Real dmin = std::numeric_limits<Real>::max();   // Initialize to "infinity"
    Size next = npoints;                            // return npoints when there is no next

//    for (Size i = 0; i < nnbs; i++)
    for (Size i = 0; i < n_nbs; i++)
    {
//        const Size     n     = points.nbs[c*nnbs+i];
        const Size     n     = points.neighbors[cum_n_nbs+i];
        const Vector3D R     = points.position[n] - points.position[o];
        const Real     Z_new = R.dot(rays.direction[r]);

        if (Z_new > Z)
        {
            const Real distance_from_ray2 = R.dot(R) - Z_new*Z_new;

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
template <>
inline Real Geometry :: get_shift <CoMoving> (
    const Size  o,
    const Size  r,
    const Size  crt                          ) const
{
    return Real(1.0) - (points.velocity[crt] - points.velocity[o]).dot(rays.direction[r]);
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline Real Geometry :: get_shift <Rest> (
    const Size  o,
    const Size  r,
    const Size  crt                      ) const
{
    Size r_correct = r;

    if (r >= nrays/2)
    {
        r_correct = rays.antipod[r];
    }

    return Real(1.0) - points.velocity[crt].dot(rays.direction[r_correct]);
}


accel inline Size Geometry :: get_n_interpl (
    const Real shift_crt,
    const Real shift_nxt,
    const Real dshift_max                   ) const
{
    const Real dshift_abs = fabs (shift_nxt - shift_crt);

    if (dshift_abs > dshift_max) {return dshift_abs/dshift_max + Real(1);}
    else                         {return 1;                              }
}



template <Frame frame>
accel inline Size Geometry :: get_ray_length (
    const Size o,
    const Size r,
    const Real dshift_max                    ) const
{
    Size  l = 0;     // ray length
    Real  Z = 0.0;   // distance from origin (o)
    Real dZ = 0.0;   // last increment in Z

    Size nxt = get_next (o, r, o, Z, dZ);

    if (nxt != get_npoints()) // if we are not going out of mesh
    {
        Size       crt = o;
        Real shift_crt = get_shift <frame> (o, r, crt);
        Real shift_nxt = get_shift <frame> (o, r, nxt);

        l += get_n_interpl (shift_crt, shift_nxt, dshift_max);

        while (boundary.point2boundary[nxt] == get_npoints()) // while nxt not on boundary
        {
                  crt =       nxt;
            shift_crt = shift_nxt;
                  nxt = get_next (o, r, nxt, Z, dZ);
//            if (nxt == get_npoints()) {printf ("ERROR (nxt < 0): o = %ld, crt = %ld, ray = %ld\n", o, crt, r); throw "";};
            shift_nxt = get_shift <frame> (o, r, nxt);

            l += get_n_interpl (shift_crt, shift_nxt, dshift_max);
        }
    }

    return l;
}



//__global__ inline void kernel (Geometry* geometry)
//{
//    printf ("Hi, from the GPU!\n");
//
//    geometry->points.position[0] += geometry->points.position[1];
//};


//inline void Geometry :: test ()
//{
//    points.position.vec[0].print();
//    points.position.vec[1].print();
//
//    Geometry* geometry_copy = (Geometry*) pc::accelerator::malloc(sizeof(Geometry));
//    pc::accelerator::memcpy_to_accelerator (geometry_copy, this, sizeof(Geometry));
//
//
//    kernel<<<1,1>>> (geometry_copy);
//
//    points.position.vec[0].print();
//    points.position.copy_ptr_to_vec ();
//    points.position.vec[0].print();
//}



inline Size1 Geometry :: get_ray_lengths ()
{
    const Size hnrays  = rays  .get_nrays  ()/2;
    const Size npoints = points.get_npoints();

    for (Size rr = 0; rr < hnrays; rr++)
    {
        const Size ar = rays.antipod[rr];

        cout << "rr = " << rr << endl;

        threaded_for (o, npoints,
        {
            const Real dshift_max = 1.0e+99;

            lengths.vec[npoints*rr+o] =   get_ray_length <CoMoving> (o, rr, dshift_max)
                                        + get_ray_length <CoMoving> (o, ar, dshift_max);
        })
    }

    return lengths.vec;
}


inline Size1 Geometry :: get_ray_lengths_gpu (
    const Size nblocks,
    const Size nthreads                      )
{
    const Size hnrays  = rays  .get_nrays  ()/2;
    const Size npoints = points.get_npoints();

    Geometry* this_cpy = (Geometry*) pc::accelerator::malloc(sizeof(Geometry));
    pc::accelerator::memcpy_to_accelerator (this_cpy, this, sizeof(Geometry));

    for (Size rr = 0; rr < hnrays; rr++)
    {
        const Size ar = rays.antipod.vec[rr];

        cout << "rr = " << rr << endl;

        accelerated_for (o, 10000, nblocks, nthreads,
        {
            const Real dshift_max = 1.0e+99;

            this_cpy->lengths[npoints*rr+o] =   this_cpy->get_ray_length <CoMoving> (o, rr, dshift_max)
                                              + this_cpy->get_ray_length <CoMoving> (o, ar, dshift_max);
        })
    }

    lengths.copy_ptr_to_vec ();

    return lengths.vec;
}
