#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "frequencies/frequencies.hpp"
#include "scattering/scattering.hpp"


///  Radiation: data structure for the radiation field
//////////////////////////////////////////////////////
struct Radiation
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model

    Frequencies frequencies;                  ///< data structure containing frequency bins

    //Data stuctures pertaining to line radiative transfer computations
    Tensor<Real> I;                           ///< intensity (r, p, f)
    Tensor<Real> u;                           ///< intensity (r, p, f)
    Tensor<Real> v;                           ///< intensity (r, p, f)

    Matrix<Real> J;                           ///< (angular) mean intensity (p, f)

    // vector<Matrix<Real>> U;         ///< U scattered intensity   (r, index(p,f))
    // vector<Matrix<Real>> V;         ///< V scattered intensity   (r, index(p,f))


    Radiation (std::shared_ptr<Parameters> params)
    : parameters  (params)
    , frequencies (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;

    inline Size index (const Size p, const Size f) const;
    inline Size index (const Size p, const Size f, const Size m) const;

    // inline Real get_U (const Size R, const Size p, const Size f) const;
    // inline Real get_V (const Size R, const Size p, const Size f) const;

    void initialize_J ();
    void MPI_reduce_J ();
    void calc_U_and_V ();
};


#include "radiation.tpp"
