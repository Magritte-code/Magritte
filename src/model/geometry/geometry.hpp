#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
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

    Vector <Size> lengths;

    void read  (const Io& io);
    void write (const Io& io) const;

    accel inline void get_next (
        const Size  o,
        const Size  r,
        const Size  crt,
              Size& nxt,
              Real& Z,
              Real& dZ,
              Real& shift ) const;

    accel inline Size get_next (
        const Size  o,
        const Size  r,
        const Size  crt,
              Real& Z,
              Real& dZ ) const;

    template <Frame frame>
    accel inline Real get_shift (
        const Size  o,
        const Size  r,
        const Size  crt ) const;

    accel inline Size get_n_interpl (
        const Real shift_crt,
        const Real shift_nxt,
        const Real dshift_max ) const;

    template <Frame frame>
    accel inline Size get_ray_length (
        const Size o,
        const Size r,
        const Real dshift_max ) const;

    inline Size1 get_ray_lengths     ();
    inline Size1 get_ray_lengths_gpu (const Size nblocks, const Size nthreads);

    inline void test ();

    inline bool valid_point     (const Size p) const;
    inline bool not_on_boundary (const Size p) const;
};


#include "geometry.tpp"
