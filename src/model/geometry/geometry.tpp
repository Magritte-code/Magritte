#include <cmath>
#include <limits>

///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      origin : position where the ray originates
///    @param[in]      raydir : normalized vector of the direction of the ray
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
template <>
accel inline Size Geometry ::get_next_general_geometry<Defaulttracer>(
    const Vector3D& origin, const Vector3D& raydir, const Size c, double& Z, double& dZ) const {
    const Size n_nbs     = points.n_neighbors[c];
    const Size cum_n_nbs = points.cum_n_neighbors[c];

    double dmin = std::numeric_limits<Real>::max(); // Initialize to "infinity"
    Size next   = parameters->npoints();            // return npoints when there is no next

    for (Size i = 0; i < n_nbs; i++) {
        const Size n       = points.neighbors[cum_n_nbs + i];
        const Vector3D R   = points.position[n] - origin;
        const double Z_new = R.dot(raydir);

        if (Z_new > Z) {
            const double distance_from_ray2 = R.dot(R) - Z_new * Z_new;
            // std::cout<<"n: "<<n<<" Z_new: "<<Z_new<<"
            // dist2:"<<distance_from_ray2<<std::endl;

            if (distance_from_ray2 < dmin) {
                dmin = distance_from_ray2;
                next = n;
                dZ   = Z_new - Z; // such that dZ > 0.0
            }
        }
    }

    // Update distance along ray
    Z += dZ;

    return next;
}

///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      origin : position where the ray originates
///    @param[in]      raydir : normalized vector of the direction of the ray
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
template <>
accel inline Size Geometry ::get_next_general_geometry<Imagetracer>(
    const Vector3D& origin, const Vector3D& raydir, const Size c, double& Z, double& dZ) const {
    const Size n_nbs     = points.n_neighbors[c];
    const Size cum_n_nbs = points.cum_n_neighbors[c];

    double dmin               = std::numeric_limits<Real>::max(); // Initialize to "infinity"
    Size next                 = parameters->npoints();            // return npoints when there is no next
    double maxdist_neighbors2 = 0.0;

    for (Size i = 0; i < n_nbs; i++) {
        const Size n       = points.neighbors[cum_n_nbs + i];
        const Vector3D R   = points.position[n] - origin;
        const double Z_new = R.dot(raydir);

        if (Z_new > Z) {
            const double distance_from_ray2 = R.dot(R) - Z_new * Z_new;

            if (distance_from_ray2 < dmin) {
                dmin = distance_from_ray2;
                next = n;
                dZ   = Z_new - Z; // such that dZ > 0.0
            }
        }
        const Vector3D diff_nxt_crt = points.position[n] - points.position[c];
        const double dist_neighbor2 = diff_nxt_crt.dot(diff_nxt_crt);
        if (dist_neighbor2 > maxdist_neighbors2) {
            maxdist_neighbors2 = dist_neighbor2;
        }
    }

    const Vector3D R_c                   = points.position[c] - origin;
    const double Z_c                     = R_c.dot(raydir);
    const double distance_curr_from_ray2 = R_c.dot(R_c) - Z_c * Z_c;

    // Precaution against tracing stuff along the boundary, stopping when we move
    // farther from the ray than the distance traveled using a safety factor of 2;
    // TO DO: actually implement some manner of 'surface unit vector' to determine
    // when to stop this ray; then this entire condition might be replaced by this
    // (maybe also merge with 'regular' get_next_general_geometry) we also check
    // whether the distance from the ray increases (as irregularly placed boundary
    // points might result in stopping too early)
    if ((!not_on_boundary(c)) && (2.0 * maxdist_neighbors2 < dmin) && (distance_curr_from_ray2 < dmin)) {
        return parameters->npoints();
    }

    // Update distance along ray
    Z += dZ;

    return next;
}

///  Getter for the number of the next cell on ray and its distance along ray
///  when assuming spherical symmetry and such that the positions are in
///  ascending order!
///    @param[in]      origin : origin from which the ray originates
///    @param[in]      raydir : direction of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
inline Size Geometry ::get_next_spherical_symmetry(
    const Vector3D& origin, const Vector3D& raydir, const Size c, double& Z, double& dZ) const {
    Size next;

    // Even though no default ray in spherical symmetry has a z-component, any
    // custom ray (for eg imaging) might have a z-component
    //  R sin is a bit difficult to determine directly, thus we use the trivial
    //  angle identity to determine Rsin^2
    const double Rcos        = origin.dot(raydir);
    const double R2          = origin.squaredNorm();
    const double Rsin2       = R2 - Rcos * Rcos;
    const double Rcos_plus_Z = Rcos + Z;

    if (Z < -Rcos) {
        if (c <= 0) {
            return parameters->npoints();
        }

        if (points.position[c - 1].squaredNorm() >= Rsin2) {
            next = c - 1;
            dZ   = -sqrt(points.position[next].squaredNorm() - Rsin2) - Rcos_plus_Z;
        } else {
            next = c;
            dZ   = -2.0 * Rcos_plus_Z;
        }
    } else {
        if (c >= parameters->npoints() - 1) {
            return parameters->npoints();
        }

        next = c + 1;
        dZ   = +sqrt(points.position[next].squaredNorm() - Rsin2) - Rcos_plus_Z;
    }

    // Update distance along ray
    Z += dZ;

    // cout << "o = " << o << "    r = " << r << "   c = " << c << "   next = " <<
    // next << "dZ = " << dZ << endl;

    return next;
}

///  Computes the distance of the ray origin to the specified boundary point in
///  any geometry. Able to both handle 3D general geometry and 1D spherical
///  symmetry.
///    @param[in]      origin : origin from which the ray originates
///    @param[in]      raydir : direction of the ray along which we are looking
///    @param[in]      bdy : point index of boundary point
inline double Geometry ::get_distance_origin_to_boundary(
    const Vector3D& origin, const Vector3D& raydir, const Size bdy) const {
    double Z = raydir.dot(points.position[bdy] - origin); // for a general 3D geometry, computing the distance is simple

    if (parameters->spherical_symmetry()) { // in spherical symmetry, the distance
                                            // needs to be computed differently;
                                            // TODO? put in geometry instead
        const double Rcos       = origin.dot(raydir);
        const double R2         = origin.squaredNorm();
        const double Rsin2      = R2 - Rcos * Rcos;
        const double horz_dist2 = points.position[bdy].squaredNorm() - Rsin2;
        // If a ray falls outside of the spherically symmetric model, we just return
        // the distance until the closest point to the model
        if (horz_dist2 > 0) {
            Z = -sqrt(horz_dist2) - Rcos; // roughly copied from geometry::get_next_spherical_symmetry
        } else {
            // Just put the distance in the middle, as the ray traced is outside of
            // the domain; this corresponds to the large 'else' clause in
            // get_next_spherical_symmetry. Then no next point will be found. We
            // expect this, as we are tracing a ray outside the model.
            Z = -Rcos;
        }
    }
    return Z;
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry ::get_shift_general_geometry<CoMoving>(
    const Vector3D& origin_velocity, const Vector3D& raydir, const Size crt) const {
    return 1.0 + (points.velocity[crt] - origin_velocity).dot(raydir);
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry ::get_shift_general_geometry<Rest>(
    const Vector3D& origin_velocity, const Vector3D& raydir, const Size crt) const {
    return 1.0 + points.velocity[crt].dot(raydir);
}

accel inline Size Geometry ::get_n_interpl(
    const double shift_crt, const double shift_nxt, const double dshift_max) const {
    const double dshift_abs = fabs(shift_nxt - shift_crt);

    if (dshift_abs > dshift_max) {
        return dshift_abs / dshift_max + 1;
    } else {
        return 1;
    }
}

template <Frame frame>
accel inline Size Geometry ::get_ray_length(const Size o, const Size r, const double dshift_max) const {
    Size l    = 0;   // ray length
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = get_next(o, r, o, Z, dZ);

    if (valid_point(nxt)) {
        Size crt         = o;
        double shift_crt = get_shift<frame>(o, r, crt, 0.0);
        double shift_nxt = get_shift<frame>(o, r, nxt, Z);

        l += 1; // no interpolation means only a single point added to the ray each
                // time

        while (not_on_boundary(nxt)) {
            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = get_next(o, r, nxt, Z, dZ);
            shift_nxt = get_shift<frame>(o, r, nxt, Z);

            l += 1; // no interpolation means only a single point added to the ray
                    // each time

            if (!valid_point(nxt)) {
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
inline bool Geometry ::valid_point(const Size p) const { return (p < parameters->npoints()); }

///  Check whether a point is not on the boundary
///    @param[in] p : point index
///    @returns true if p is not on the boundary
/////////////////////////////////////////////////
inline bool Geometry ::not_on_boundary(const Size p) const {
    return (boundary.point2boundary[p] == parameters->npoints());
}

///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      origin : position from which the ray originates
///    @param[in]      raydir : direction of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
template <Tracer tracer>
accel inline Size Geometry ::get_next(
    const Vector3D& origin, const Vector3D& raydir, const Size crt, double& Z, double& dZ) const {
    Size next;

    if (parameters->spherical_symmetry()) {
        next = get_next_spherical_symmetry(origin, raydir, crt, Z, dZ);
    } else {
        next = get_next_general_geometry<tracer>(origin, raydir, crt, Z, dZ);
    }

    return next;
}

///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///  Simple wrapper for compatibility with the previous internal api
///    @param[in]      o : number of cell from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
accel inline Size Geometry ::get_next(const Size o, const Size r, const Size crt, double& Z, double& dZ) const {
    const Vector3D origin = points.position[o];
    const Vector3D raydir = rays.direction[r];

    return get_next<Defaulttracer>(origin, raydir, crt, Z, dZ);
}

///  Getter for the number of the next cell on ray and its distance along ray in
///  the general case without any further assumptions
///    @param[in]      origin : position from which the ray originates
///    @param[in]      raydir : direction of the ray along which we are looking
///    @param[in]      c : number of the cell put last on the ray
///    @param[in/out]  Z : reference to the current distance along the ray
///    @param[out]    dZ : reference to the distance increment to the next ray
///    @return number of the next cell on the ray after the current cell
///////////////////////////////////////////////////////////////////////////////////
accel inline void Geometry ::get_next(
    const Size o, const Size r, const Size crt, Size& nxt, double& Z, double& dZ, double& shift) const {
    nxt   = get_next(o, r, crt, Z, dZ);
    shift = get_shift<CoMoving>(o, r, nxt, Z);
}

///  Getter for the number of the next closer boundary point on the ray and its
///  distance along ray in the general case without any further assumptions. Can
///  handle rays not originating from inside the model. Remark: we assume the
///  outer boundary to be convex (as usual in Magritte) in order to find a
///  consistent result when applying this multiple times to find the closest
///  boundary point to the given ray
///    @param[in] origin : coordinate of the origin point from which the ray
///    originates
///    @param[in] raydir : direction vector of the ray along which we are
///    looking
///    @param[in] crt : number of the cell put last on the ray
///    @param[in/out] Z : reference to the current distance along the ray
///    @param[out] dZ : reference to the distance increment to the next ray
///    @return number of the next closer boundary point on the ray
//////////////////////////////////////////////////////////////////
accel inline Size Geometry ::get_boundary_point_closer_to_custom_ray(
    const Vector3D& origin, const Vector3D& raydir, const Size crt) const
// double& Z,
// double& dZ      ) const
{
    const Size n_nbs     = points.n_neighbors[crt];
    const Size cum_n_nbs = points.cum_n_neighbors[crt];

    const Vector3D R_curr = points.position[crt] - origin;
    const double Z_curr   = R_curr.dot(raydir);

    double dmin2 = R_curr.dot(R_curr) - Z_curr * Z_curr; // Initialize to current distance
    Size next    = crt;                                  // return current point when no better boundary point is found

    for (Size i = 0; i < n_nbs; i++) {
        const Size n = points.neighbors[cum_n_nbs + i];
        // TODO: Lazy implementation: define instead a neighbors structure using
        // only boundary points (now we are looping over way too many points)
        if (!not_on_boundary(n)) {
            const Vector3D R_new            = points.position[n] - origin;
            const double Z_new              = R_new.dot(raydir);
            const double distance_from_ray2 = R_new.dot(R_new) - Z_new * Z_new;

            if (distance_from_ray2 < dmin2) {
                dmin2 = distance_from_ray2;
                next  = n;
            }
        }
    }

    return next;
}

///  Compute the closest boundary point in the custom ray direction
///    @param[in] raydir : ray direction vector
///////////////////////////////////////////////////////////////////
accel inline Size Geometry ::get_closest_bdy_point_in_custom_raydir(const Vector3D& raydir) const {
    // in spherical symmetry, the furthest point is always the last point
    if (parameters->spherical_symmetry()) {
        return parameters->npoints() - 1;
    }
    // first try out first boundary point
    Size closest_bdy_point = boundary.point2boundary[0]; // first best guess is the first boundary point
    double projected_dmin  = raydir.dot(points.position[closest_bdy_point]); // and the corresponding
                                                                             // projected distance
    for (Size bdy_index = 1; bdy_index < parameters->nboundary(); bdy_index++) {
        const Size bdy_point_index      = boundary.boundary2point[bdy_index];
        const double projected_distance = raydir.dot(points.position[bdy_point_index]);

        if (projected_distance < projected_dmin) {
            closest_bdy_point = bdy_point_index;
            projected_dmin    = projected_distance;
        }
    }

    return closest_bdy_point;
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry ::get_shift_spherical_symmetry<CoMoving>(const Vector3D& origin, const Vector3D& origin_velocity,
    const Vector3D& raydir, const Size c, const double Z) const {
    // At the 0 point, we cannot have any velocity (otherwise not spherically
    // symmetric)
    if (origin.squaredNorm() == 0.0) {
        return 1.0;
    }

    if (points.position[c].x() == 0.0) {
        return 1.0 - origin_velocity.dot(raydir);
    }

    const double Rcos_plus_Z = origin.dot(raydir) + Z;

    return 1.0 + (points.velocity[c].x() * Rcos_plus_Z / points.position[c].x() - origin_velocity.dot(raydir));
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <>
inline double Geometry ::get_shift_spherical_symmetry<Rest>(const Vector3D& origin, const Vector3D& origin_velocity,
    const Vector3D& raydir, const Size c, const double Z) const {
    if (points.position[c].x() == 0.0) {
        return 1.0;
    }

    const double Rcos_plus_Z = origin.dot(raydir) + Z;

    return 1.0 + points.velocity[c].x() * Rcos_plus_Z / points.position[c].x();
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin Wrapper for compatibility with old api
///    @param[in] o   : number of cell from which the ray originates
///    @param[in] r   : number of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <Frame frame>
inline double Geometry ::get_shift(const Size o, const Size r, const Size c, const double Z) const {
    Vector3D raydir                = rays.direction[r];
    const Vector3D origin          = points.position[o];
    const Vector3D origin_velocity = points.velocity[o];
    // Due to the old imager implementation, the computed shift might need to be
    // reversed
    const bool reverse = ((frame == Rest) && (r >= parameters->hnrays()));

    return get_shift<frame>(origin, origin_velocity, raydir, c, Z, reverse);
}

///  Getter for the doppler shift along the ray between the current cell and the
///  origin
///    @param[in] origin : position from which the ray originates
///    @param[in] raydir : direction of the ray along which we are looking
///    @param[in] crt : number of the cell for which we want the velocity
///    @param[in] Z   : reference to the current distance along the ray
///    @param[in] reverse: whether to reverse the computed doppler shift (for
///    imaging)
///    @return doppler shift along the ray between the current cell and the
///    origin
///////////////////////////////////////////////////////////////////////////////////////
template <Frame frame>
inline double Geometry ::get_shift(const Vector3D& origin, const Vector3D& origin_velocity, const Vector3D& raydir,
    const Size c, const double Z, const bool reverse) const {
    double shift = 1.0;
    if (parameters->spherical_symmetry()) {
        shift = get_shift_spherical_symmetry<frame>(origin, origin_velocity, raydir, c, Z);
    } else {
        shift = get_shift_general_geometry<frame>(origin_velocity, raydir, c);
    }
    if (reverse) {
        return 2.0 - shift;
    } else {
        return shift;
    }
}
