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

    Tensor<Real> I;         ///< intensity (r, p, f)
    Matrix<Real> J;         ///< (angular) mean intensity (p, f)

    Real1 print(Size r, Size o)
    {
        for (Size f = 0; f < parameters.nfreqs(); f++)
        {
            printf("I(f=%ld) = %le\n", f, I(r,o,f));

            return I.vec;
        }
    }

    // vector<Matrix<Real>> u;         ///< u intensity             (r, index(p,f))
    // vector<Matrix<Real>> v;         ///< v intensity             (r, index(p,f))

    // vector<Matrix<Real>> U;         ///< U scattered intensity   (r, index(p,f))
    // vector<Matrix<Real>> V;         ///< V scattered intensity   (r, index(p,f))

    // Real1 J;         ///< (angular) mean intensity (index(p,f))
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

    void initialize_J ();
    void MPI_reduce_J ();
    void calc_U_and_V ();
};


#include "radiation.tpp"
