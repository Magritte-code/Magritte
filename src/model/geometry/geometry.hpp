#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "tools/constants.hpp"
#include "points/points.hpp"
#include "rays/rays.hpp"
#include "boundary/boundary.hpp"

///  Frame of reference used in geometry computations
/////////////////////////////////////////////////////
enum Frame {CoMoving, Rest};


///  Data structure for geometry
////////////////////////////////
struct Geometry
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    Points   points;                          ///< data structure containing point data
    Rays     rays;                            ///< data structure containing ray (direction) data
    Boundary boundary;                        ///< data structure containing boundary data

    Matrix<Size> lengths;
    Size         lengths_max;


    Geometry (std::shared_ptr<Parameters> params)
    : parameters (params)
    , points     (params)
    , rays       (params)
    , boundary   (params) {};


    void read  (const Io& io);
    void write (const Io& io) const;

    accel inline void get_next (
        const Size    o,
        const Size    r,
        const Size    crt,
              Size&   nxt,
              double& Z,
              double& dZ,
              double& shift ) const;

    accel inline Size get_next (
        const Size    o,
        const Size    r,
        const Size    crt,
              double& Z,
              double& dZ  ) const;

    // accel inline Size get_next_custom_origin_raydir (
    //     const Vector3D origin,
    //     const Vector3D raydir,
    //     const Size      crt,
    //           double&   Z,
    //           double&   dZ    ) const;

    accel inline Size get_next_general_geometry (
        const Size    o,
        const Size    r,
        const Size    crt,
              double& Z,
              double& dZ  ) const;

    accel inline Size get_next_general_geometry_custom_origin_raydir (
        const Vector3D origin,
        const Vector3D raydir,
        const Size     crt,
              double&  Z,
              double&  dZ     ) const;

    accel inline Size get_next_spherical_symmetry (
        const Size    o,
        const Size    r,
        const Size    crt,
              double& Z,
              double& dZ  ) const;

    // accel inline Size get_next_spherical_symmetry_custom_origin_raydir (
    //     const Vector3D origin,
    //     const Vector3D raydir,
    //     const Size     crt,
    //           double&  Z,
    //           double&  dZ     ) const;

    accel inline Size get_boundary_point_closer_to_custom_ray (
        const Vector3D origin,
        const Vector3D raydir,
        const Size     crt    ) const;

    accel inline Size get_closest_bdy_point_in_custom_raydir (
        const Vector3D raydir) const;

    template <Frame frame>
    accel inline double get_shift (
        const Size   o,
        const Size   r,
        const Size   crt,
        const double Z   ) const;

    accel inline double get_shift_general_geometry_custom_origin_raydir (
        const Vector3D raydir,
        const Size     crt) const;

    template <Frame frame>
    accel inline double get_shift_general_geometry (
        const Size o,
        const Size r,
        const Size crt ) const;

    template <Frame frame>
    accel inline double get_shift_spherical_symmetry (
        const Size   o,
        const Size   r,
        const Size   crt,
        const double Z   ) const;

    accel inline Size get_n_interpl (
        const double shift_crt,
        const double shift_nxt,
        const double dshift_max ) const;

    template <Frame frame>
    accel inline Size get_ray_length (
        const Size   o,
        const Size   r,
        const double dshift_max ) const;

    inline bool valid_point     (const Size p) const;
    inline bool not_on_boundary (const Size p) const;
};


#include "geometry.tpp"
