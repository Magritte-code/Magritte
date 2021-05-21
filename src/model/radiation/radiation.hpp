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
    Tensor<Real> u;         ///< intensity (r, p, f)
    // Tensor<Real> v;         ///< intensity (r, p, f)

    Matrix<Real> J;         ///< (angular) mean intensity (p, f)

    // vector<Matrix<Real>> U;         ///< U scattered intensity   (r, index(p,f))
    // vector<Matrix<Real>> V;         ///< V scattered intensity   (r, index(p,f))

    void read  (const Io& io);
    void write (const Io& io) const;

    inline Size index (const Size p, const Size f) const;
    inline Size index (const Size p, const Size f, const Size m) const;

    void initialize_J ();
    void MPI_reduce_J ();
    void calc_U_and_V ();
};


#include "radiation.tpp"
