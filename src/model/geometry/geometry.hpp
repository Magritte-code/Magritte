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
    Parameters parameters;
    Points     points;
    Rays       rays;
    Boundary   boundary;

    Matrix<Size> lengths;
    Size         lengths_max;

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

    accel inline Size get_next_general_geometry (
        const Size    o,
        const Size    r,
        const Size    crt,
              double& Z,
              double& dZ  ) const;

    accel inline Size get_next_spherical_symmetry (
        const Size    o,
        const Size    r,
        const Size    crt,
              double& Z,
              double& dZ  ) const;

    template <Frame frame>
    accel inline double get_shift (
        const Size   o,
        const Size   r,
        const Size   crt,
        const double Z   ) const;

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

    // template <Frame frame>
    // inline void get_ray_lengths (const double dshift_max);

    inline bool valid_point     (const Size p) const;
    inline bool not_on_boundary (const Size p) const;
};


#include "geometry.tpp"
