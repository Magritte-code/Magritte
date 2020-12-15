#include <limits>
#include <cmath>
#include <algorithm>    // std::max
#include <set>


///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      o : number of cell from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
accel inline Size Geometry :: get_next_general_geometry (
    const Size    o,
    const Size    r,
    const Size    c,
          double& Z,
          double& dZ                   ) const
{
    Size n_nbs = points.multiscale.get_nb_neighbors(c);//    n_neighbors[c];
//    const Size cum_n_nbs = points.cum_n_neighbors[c];

    double dmin = std::numeric_limits<Real>::max();   // Initialize to "infinity"
    Size   next = parameters.npoints();               // return npoints when there is no next

//    for (Size i = 0; i < nnbs; i++)

    //TODO: update to use set instead of vector
    std::set<Size> temp_neighbors=points.multiscale.get_neighbors(c);
    // temp_neighbors.insert(std::end(temp_neighbors), std::begin(points.multiscale.get_neighbors(c).begin()), std::end(points.multiscale.get_neighbors(c).end()));
    // Vector<Size> temp_neighbors(temp_vector);
    for (Size n:temp_neighbors)
    // for (Size i = 0; i < n_nbs; i++)
    {
//        const Size     n     = points.nbs[c*nnbs+i];
        // const Size     n     = temp_neighbors[i];//points.neighbors[cum_n_nbs+i];
        const Vector3D R     = points.position[n] - points.position[o];
        const double   Z_new = R.dot(rays.direction[r]);

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


///  Getter for the number of the next cell on ray and its distance along ray when
///  assuming spherical symmetry and such that the positions are in ascending order!
///    @param[in]      o : number of cell from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
inline Size Geometry :: get_next_spherical_symmetry (
    const Size  o,
    const Size  r,
    const Size  c,
    double     &Z,
    double     &dZ                                  ) const
{
    Size next;

    const double Rsin = points.position[o].x() * rays.direction[r].y();
    const double Rcos = points.position[o].x() * rays.direction[r].x();

    const double Rsin2       = Rsin * Rsin;
    const double Rcos_plus_Z = Rcos + Z;

    if (Z < -Rcos)
    {
        if (c <= 0)
        {
            return parameters.npoints();
        }

        if (points.position[c-1].squaredNorm() >= Rsin2)
        {
            next = c - 1;
            Size curr_coars_lvl=points.multiscale.get_curr_coars_lvl();
            while(!(points.multiscale.get_mask(curr_coars_lvl))[next])
            {
              // std::cout<<"next: "<<next<<std::endl;
              next=next-1;
            }
            dZ   = -sqrt(points.position[next].squaredNorm() - Rsin2) - Rcos_plus_Z;
        }
        else
        {
            next = c;
            dZ   = - 2.0 * Rcos_plus_Z;
        }
    }
    else
    {
        if (c >= parameters.npoints()-1)
        {
            return parameters.npoints();
        }

        next = c + 1;
        Size curr_coars_lvl=points.multiscale.get_curr_coars_lvl();
        while(!(points.multiscale.get_mask(curr_coars_lvl))[next])
        {
          // std::cout<<"next: "<<next<<std::endl;
          next=next+1;
        }
        dZ   = +sqrt(points.position[next].squaredNorm() - Rsin2) - Rcos_plus_Z;
    }

    // Update distance along ray
    Z += dZ;

    // cout << "o = " << o << "    r = " << r << "   c = " << c << "   next = " << next << "dZ = " << dZ << endl;

    return next;
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry :: get_shift_general_geometry <CoMoving> (
    const Size  o,
    const Size  r,
    const Size  crt ) const
{
    return 1.0 - (points.velocity[crt] - points.velocity[o]).dot(rays.direction[r]);
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry :: get_shift_general_geometry <Rest> (
    const Size  o,
    const Size  r,
    const Size  crt ) const
{
    Size r_correct = r;

    if (r >= parameters.hnrays()) // assumes ray indices and antipodes are on opposite sites of hnrays
    {
        r_correct = rays.antipod[r];
    }

    return 1.0 - points.velocity[crt].dot(rays.direction[r_correct]);
}


accel inline Size Geometry :: get_n_interpl (
    const double shift_crt,
    const double shift_nxt,
    const double dshift_max                 ) const
{
    const double dshift_abs = fabs (shift_nxt - shift_crt);

    if (dshift_abs > dshift_max) {return dshift_abs/dshift_max + 1;}
    else                         {return 1;                        }
}



template <Frame frame>
accel inline Size Geometry :: get_ray_length (
    const Size   o,
    const Size   r,
    const double dshift_max                  ) const
{
    Size    l = 0;     // ray length
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = get_next (o, r, o, Z, dZ);

    if (valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = get_shift <frame> (o, r, crt, Z);
        double shift_nxt = get_shift <frame> (o, r, nxt, Z);

        l += get_n_interpl (shift_crt, shift_nxt, dshift_max);

        while (not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = get_next          (o, r, nxt, Z, dZ);
            shift_nxt = get_shift <frame> (o, r, nxt, Z    );

            l += get_n_interpl (shift_crt, shift_nxt, dshift_max);
        }
    }

    return l;
}


// inline Size1 Geometry :: get_ray_lengths ()
// {
//     for (Size rr = 0; rr < parameters.hnrays(); rr++)
//     {
//         const Size ar = rays.antipod[rr];
//
//         cout << "rr = " << rr << endl;
//
//         threaded_for (o, parameters.npoints(),
//         {
//             const Real dshift_max = 1.0e+99;
//
//             lengths.vec[parameters.npoints()*rr+o] =
//                 get_ray_length <CoMoving> (o, rr, dshift_max)
//               + get_ray_length <CoMoving> (o, ar, dshift_max);
//         })
//     }
//
//     return lengths.vec;
// }


// template <Frame frame>
// inline void Geometry :: get_ray_lengths (const Real dshift_max)
// {
//     for (Size rr = 0; rr < parameters.hnrays(); rr++)
//     {
//         const Size ar = rays.antipod[rr];
//
//         accelerated_for (o, parameters.npoints(), nblocks, nthreads,
//         {
//             lengths(rr,o) =  get_ray_length <frame> (o, rr, dshift_max)
//                            + get_ray_length <frame> (o, ar, dshift_max);
//         })
//
//         pc::accelerator::synchronize();
//     }
//
//     lengths.copy_ptr_to_vec();
// }


///  Check whether a point index is valid
///    @param[in] p : point index
///    @returns true if p is a valid index
//////////////////////////////////////////
inline bool Geometry :: valid_point (const Size p) const
{
    return (p < parameters.npoints());
}


///  Check whether a point is not on the boundary
///    @param[in] p : point index
///    @returns true if p is not on the boundary
/////////////////////////////////////////////////
inline bool Geometry :: not_on_boundary (const Size p) const
{
    return (boundary.point2boundary[p] == parameters.npoints());
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
    const Size    o,
    const Size    r,
    const Size    crt,
          double& Z,
          double& dZ                   ) const
{
    Size next;

    if (parameters.spherical_symmetry())
    {
        next = get_next_spherical_symmetry (o, r, crt, Z, dZ);
    }
    else
    {
        next = get_next_general_geometry   (o, r, crt, Z, dZ);
    }

    //if (!valid_point (next))
    //{
    //    printf ("ERROR (next is not valid): o = %d, crt = %d, ray = %d\n", o, crt, r);
    //}

    return next;
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
accel inline void Geometry :: get_next (
    const Size    o,
    const Size    r,
    const Size    crt,
          Size&   nxt,
          double& Z,
          double& dZ,
          double& shift ) const
{
    nxt   = get_next             (o, r, crt, Z, dZ);
    shift = get_shift <CoMoving> (o, r, nxt, Z    );
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template<>
inline double Geometry :: get_shift_spherical_symmetry <CoMoving> (
    const Size   o,
    const Size   r,
    const Size   c,
    const double Z ) const
{
    if (points.position[o].x() == 0.0)
    {
        return 1.0;
    }

    if (points.position[c].x() == 0.0)
    {
        return 1.0 + points.velocity[o].x() * rays.direction[r].x();
    }

    const double Rcos_plus_Z = points.position[o].x() * rays.direction[r].x() + Z;

    return 1.0 - (  points.velocity[c].x() * Rcos_plus_Z / points.position[c].x()
                  - points.velocity[o].x() * rays.direction[r].x()               );
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template<>
inline double Geometry :: get_shift_spherical_symmetry <Rest> (
    const Size   o,
    const Size   r,
    const Size   c,
    const double Z ) const
{
    if (points.position[c].x() == 0.0)
    {
        return 1.0;
    }

    const double Rcos_plus_Z = points.position[o].x() * rays.direction[r].x() + Z;

    // double shift;

    if (r < parameters.hnrays()) // assumes ray indices and antipodes are on opposite sites of hnrays
    {
        return 1.0 - points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    }
    else
    {
        return 1.0 + points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    }

    // cout << "r = " << r << "  o = " << o << "  c = " << c << "  " << shift << endl;

    // return shift;


    // const double shift = 1.0 - points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    // return shift;

    //if (rays.direction[r].x() >= 0)
    //{
    //    const double shift = 1.0 - points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    //    return shift;
    //}
    //else
    //{
    //    const double shift = 1.0 + points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    //    return shift;
    //}
}


///  Getter for the doppler shift along the ray between the current cell and the origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the origin
///////////////////////////////////////////////////////////////////////////////////////
template<Frame frame>
inline double Geometry :: get_shift (
    const Size   o,
    const Size   r,
    const Size   c,
    const double Z ) const
{
    if (parameters.spherical_symmetry())
    {
        return get_shift_spherical_symmetry <frame> (o, r, c, Z);
    }
    else
    {
        return get_shift_general_geometry <frame> (o, r, c);
    }
}
