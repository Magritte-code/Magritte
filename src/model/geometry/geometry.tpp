#include <limits>
#include <cmath>


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
    const Size     n_nbs = points.    n_neighbors[c];
    const Size cum_n_nbs = points.cum_n_neighbors[c];

    double dmin = std::numeric_limits<Real>::max();   // Initialize to "infinity"
    Size   next = parameters->npoints();              // return npoints when there is no next

    for (Size i = 0; i < n_nbs; i++)
    {
        const Size     n     = points.neighbors[cum_n_nbs+i];
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
            return parameters->npoints();
        }

        if (points.position[c-1].squaredNorm() >= Rsin2)
        {
            next = c - 1;
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
        if (c >= parameters->npoints()-1)
        {
            return parameters->npoints();
        }

        next = c + 1;
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

    if (r >= parameters->hnrays()) // assumes ray indices and antipodes are on opposite sites of hnrays
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
        double shift_crt = get_shift <frame> (o, r, crt, 0.0);
        double shift_nxt = get_shift <frame> (o, r, nxt, Z);

        l += 1;//no interpolation means only a single point added to the ray each time

        while (not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = get_next          (o, r, nxt, Z, dZ);
            shift_nxt = get_shift <frame> (o, r, nxt, Z    );

            l += 1;//no interpolation means only a single point added to the ray each time

            if (!valid_point(nxt))
            {
                printf("ERROR: no valid neighbor o=%u, r=%u, crt=%u\n", o, r, crt);
            }
        }
    }

    return l;
}


///  Check whether a point index is valid
///    @param[in] p : point index
///    @returns true if p is a valid index
//////////////////////////////////////////
inline bool Geometry :: valid_point (const Size p) const
{
    return (p < parameters->npoints());
}


///  Check whether a point is not on the boundary
///    @param[in] p : point index
///    @returns true if p is not on the boundary
/////////////////////////////////////////////////
inline bool Geometry :: not_on_boundary (const Size p) const
{
    return (boundary.point2boundary[p] == parameters->npoints());
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

    if (parameters->spherical_symmetry())
    {
        next = get_next_spherical_symmetry (o, r, crt, Z, dZ);
    }
    else
    {
        next = get_next_general_geometry   (o, r, crt, Z, dZ);
    }

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


///  Getter for the number of the next closer boundary point on the ray and its distance along ray in
///  the general case without any further assumptions. Can handle rays not originating from inside the model.
///  Remark: we assume the outer boundary to be convex (as usual in Magritte) in order to find a consistent result when applying this multiple times to find the closest boundary point to the given ray
///    @param[in] origin : coordinate of the origin point from which the ray originates
///    @param[in] raydir : direction vector of the ray along which we are looking
///    @param[in] crt : number of the cell put last on the ray
///    @param[in/out] Z : reference to the current distance along the ray
///    @param[out] dZ : reference to the distance increment to the next ray
///    @return number of the next closer boundary point on the ray
//////////////////////////////////////////////////////////////////
accel inline Size Geometry :: get_boundary_point_closer_to_custom_ray (
    const Vector3D origin,
    const Vector3D raydir,
    const Size     crt) const
          // double& Z,
          // double& dZ      ) const
{
    const Size     n_nbs = points.    n_neighbors[crt];
    const Size cum_n_nbs = points.cum_n_neighbors[crt];

    const Vector3D R_curr = points.position[crt] - origin;
    const double   Z_curr = R_curr.dot(raydir);

    double dmin2 = R_curr.dot(R_curr) - Z_new*Z_new;   // Initialize to current distance
    Size   next = crt;              // return current point when no better boundary point is found

    for (Size i = 0; i < n_nbs; i++)
    {
        const Size     n     = points.neighbors[cum_n_nbs+i];
        //TODO: Lazy implementation: define instead a neighbors structure using only boundary points (now we are looping over way too many points)
        if (!not_on_boundary(n))
        {
            const Vector3D R_new = points.position[n] - origin;
            const double   Z_new = R.dot(raydir);
            const double distance_from_ray2 = R_new.dot(R_new) - Z_new*Z_new;

            if (distance_from_ray2 < dmin2)
            {
                dmin2 = distance_from_ray2;
                next = n;
            }
        }
    }

    return next;
}


///  Compute the closest boundary point in the custom ray direction
///    @param[in] raydir : ray direction vector
///////////////////////////////////////////////////////////////////
accel inline Size Geometry :: get_closest_bdy_point_in_custom_raydir (
    const Vector3D raydir) const
{
    //first try out first boundary point
    Size closest_bdy_point = boundary.point2boundary[0];//first best guess is the first boundary point
    double projected_dmin = raydir.dot(points.position[closest_point]); // and the corresponding projected distance
    for (Size bdy_index = 1; bdy_index<parameters->nboundary(); bdy_index++)
    {
        const Size bdy_point_index = boundary.point2boundary[bdy_index];
        const double projected_distance = raydir.dot(points.position[bdy_point_index]);

        if (projected_distance<projected_dmin)
        {
            closest_bdy_point = bdy_point_index;
        }
    }

    return closest_bdy_point;
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

    if (r < parameters->hnrays()) // assumes ray indices and antipodes are on opposite sites of hnrays
    {
        return 1.0 - points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    }
    else
    {
        return 1.0 + points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
    }
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
    if (parameters->spherical_symmetry())
    {
        return get_shift_spherical_symmetry <frame> (o, r, c, Z);
    }
    else
    {
        return get_shift_general_geometry <frame> (o, r, c);
    }
}
