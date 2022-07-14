#pragma once

#include <Eigen/SparseLU>
using Eigen::SparseLU;
#include <Eigen/SparseCore>
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::Triplet;
using Eigen::COLAMDOrdering;


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "linedata/linedata.hpp"
#include "quadrature/quadrature.hpp"
#include "lambda/lambda.hpp"


struct LineProducingSpecies
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    Linedata   linedata;                      ///< data for line producing species
    Quadrature quadrature;                    ///< data for integral over line
    Lambda     lambda;                        ///< Approximate Lambda Operator (ALO)

    Real2 Jlin;                               ///< actual mean intensity in the line
    Real2 Jeff;                               ///< effective mean intensity in the line (actual - ALO)
    Real2 Jdif;                               ///< effective mean intensity in the line (actual - ALO)

    Matrix<Real> J;                           ///< Isotropic radiation field
    Matrix<Real> J2_0;                        ///< Anisotropic radiation field tensor element 0
    Matrix<Real> J2_1_Re;                     ///< Anisotropic radiation field tensor element 1 (real part)
    Matrix<Real> J2_1_Im;                     ///< Anisotropic radiation field tensor element 1 (imaginary part)
    Matrix<Real> J2_2_Re;                     ///< Anisotropic radiation field tensor element 2 (real part)
    Matrix<Real> J2_2_Im;                     ///< Anisotropic radiation field tensor element 2 (imaginary part)

    Size3 nr_line;                            ///< frequency number corresponing to line (p,k,z)

    double relative_change_mean;              ///< mean    relative change
    double relative_change_max;               ///< maximum relative change
    double fraction_not_converged;            ///< fraction of levels that is not converged

    VectorXr population;                      ///< level population (most recent)
    Real1    population_tot;                  ///< total level population (sum over levels)

    vector<VectorXr> populations;             ///< list of populations in previous iterations
    vector<VectorXr> residuals;               ///< list of residuals in the populations


    VectorXr population_prev1;                ///< level populations 1 iteration  back
    VectorXr population_prev2;                ///< level populations 2 iterations back
    VectorXr population_prev3;                ///< level populations 3 iterations back

    SparseMatrix<Real> RT;
    SparseMatrix<Real> LambdaTest;
    SparseMatrix<Real> LambdaStar;


    // PORTAL extension
    ////////////////////////////////
    Size  nalign;
    Size1 j_lev;
    Size2 a_lev;
    Double2 J0;
    Double2 J2;
    Double5 rp;
    Double5 rm;
    Double2 tp;

    MatrixXd M;

    Double2 rho;
    ////////////////////////////////


    LineProducingSpecies (std::shared_ptr<Parameters> params)
    : parameters (params)
    , quadrature (params)
    , lambda     (params) {};

    void read  (const Io& io, const Size l);
    void write (const Io& io, const Size l) const;

    void read_populations  (const Io& io, const Size l, const string tag);
    void write_populations (const Io& io, const Size l, const string tag) const;

    inline Size index (const Size p, const Size i) const;

    inline Real get_emissivity (const Size p, const Size k) const;
    inline Real get_opacity    (const Size p, const Size k) const;

    inline void check_for_convergence (
        const Real pop_prec );

    inline void update_using_LTE (
        const Double2      &abundance,
        const Vector<Real> &temperature );

    inline void update_using_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature );

    inline void update_using_statistical_equilibrium_sparse (
        const Double2      &abundance,
        const Vector<Real> &temperature );

    inline void update_using_Ng_acceleration ();
    inline void update_using_acceleration (const Size order);

    inline void PORTAL_solve_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature );

    inline void PORTAL_solve_statistical_equilibrium_for_point (
        const Double2      &abundance,
        const Vector<Real> &temperature,
        const Size          p           );
};


#include "lineProducingSpecies.tpp"
