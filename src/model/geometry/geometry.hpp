#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "points/points.hpp"
#include "rays/rays.hpp"
#include "boundary/boundary.hpp"

///  Frame of reference used in geometry computations
/////////////////////////////////////////////////////
enum Frame {CoMoving, Rest};


///  Geometry class
///////////////////
class Geometry
{
    public:
        Points   points;
        Rays     rays;
        Boundary boundary;

        void read  (const Io& io);
        void write (const Io& io) const;

        inline void   set_npoints (const size_t n);
        inline size_t get_npoints () const;

//        inline void trace_ray (const size_t origin, const size_t ray) const;

        inline size_t get_next (
            const size_t  o,
            const size_t  r,
            const size_t  c,
                  double& Z,
                  double& dZ   ) const;

        template <Frame frame>
        inline double get_shift (
            const size_t  o,
            const size_t  r,
            const size_t  crt   ) const;


        inline size_t get_n_interpl (
            const double  shift_crt,
            const double  shift_nxt,
            const double  dshift_max ) const;

        template <Frame frame>
        inline size_t get_ray_length (
            const size_t    o,
            const size_t    r,
            const double    dshift_max ) const;

        inline Long2 get_ray_lengths () const;

    private:
        size_t npoints;
        size_t nrays;
};


#include "geometry.tpp"