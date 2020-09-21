#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "frequencies/frequencies.hpp"
#include "scattering/scattering.hpp"


///  Radiation: data structure for the radiation field
//////////////////////////////////////////////////////
struct Radiation
{
    Parameters parameters;

    Frequencies frequencies;
    Scattering  scattering;

    Real2 u;         ///< u intensity             (r, index(p,f))
    Real2 v;         ///< v intensity             (r, index(p,f))

    Real2 U;         ///< U scattered intensity   (r, index(p,f))
    Real2 V;         ///< V scattered intensity   (r, index(p,f))

    Real1 J;         ///< (angular) mean intensity (index(p,f))

    Real3 I_bdy;     ///< intensity at the boundary (r,b,f)

    void read  (const Io& io);
    void write (const Io& io) const;


    inline Size index (const Size p, const Size f) const;
    inline Size index (const Size p, const Size f, const Size m) const;

    inline Real get_U (const Size R, const Size p, const Size f) const;
    inline Real get_V (const Size R, const Size p, const Size f) const;

    inline Real get_I_bdy (const Size R, const Size p, const Size f) const;


    inline Real get_u (const Size r, const Size p, const Size f) const;
    inline Real get_v (const Size r, const Size p, const Size f) const;

    inline Real get_J (const Size p, const Size f) const;


    //inline void rescale_U_and_V (
    //    const Real &freq_scaled,
    //    const Size  R,
    //    const Size  p,
    //          Size &notch,
    //          Real &U_scaled,
    //          Real &V_scaled   ) const;

    //inline void rescale_I_bdy (
    //    const Real &freq_scaled,
    //    const Size  R,
    //    const Size  p,
    //    const Size  b,
    //          Size &notch,
    //          Real &Ibdy_scaled) const;


    void initialize_J ();
    void MPI_reduce_J ();
    void calc_U_and_V ();
};


#include "radiation.tpp"
