#pragma once

#include <Eigen/SparseLU>
using Eigen::SparseLU;
#include <Eigen/SparseCore>
using Eigen::COLAMDOrdering;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::Triplet;

#include "io/io.hpp"
#include "lambda/lambda.hpp"
#include "linedata/linedata.hpp"
#include "model/parameters/parameters.hpp"
#include "quadrature/quadrature.hpp"
#include "tools/types.hpp"

struct LineProducingSpecies {
    std::shared_ptr<Parameters> parameters; ///< data structure containing model parameters

    Linedata linedata;     ///< data for line producing species
    Quadrature quadrature; ///< data for integral over line
    Lambda lambda;         ///< Approximate Lambda Operator (ALO)

    Real2 Jlin; ///< actual mean intensity in the line
    Real2 Jeff; ///< effective mean intensity in the line (actual - ALO)
    Real2 Jdif; ///< effective mean intensity in the line (actual - ALO)

    Matrix<Real> J;       ///< Isotropic radiation field
    Matrix<Real> J2_0;    ///< Anisotropic radiation field tensor element 0
    Matrix<Real> J2_1_Re; ///< Anisotropic radiation field tensor element 1 (real part)
    Matrix<Real> J2_1_Im; ///< Anisotropic radiation field tensor element 1
                          ///< (imaginary part)
    Matrix<Real> J2_2_Re; ///< Anisotropic radiation field tensor element 2 (real part)
    Matrix<Real> J2_2_Im; ///< Anisotropic radiation field tensor element 2
                          ///< (imaginary part)

    Size3 nr_line; ///< frequency number corresponing to line (p,k,z)

    double relative_change_mean;   ///< mean    relative change
    double relative_change_max;    ///< maximum relative change
    double fraction_not_converged; ///< fraction of levels that is not converged

    // For ng-acceleration purposes, the level populations must be stored as accurately as possible

    VectorXld population; ///< level population (most recent)
    Real1 population_tot; ///< total level population (sum over levels)

    vector<VectorXld> populations; ///< list of populations in previous iterations
    vector<VectorXld> residuals;   ///< list of residuals in the populations

    VectorXld trial_population;      ///< level population generated in adaptive ng
                                     ///< acceleration trial
    VectorXld trial_population_prev; ///< level population generated in the previous adaptive ng
                                     ///< acceleration trial

    VectorXld population_prev1; ///< level populations 1 iteration  back
    VectorXld population_prev2; ///< level populations 2 iterations back
    VectorXld population_prev3; ///< level populations 3 iterations back

    SparseMatrix<long double> RT;
    SparseMatrix<long double> LambdaTest;
    SparseMatrix<long double> LambdaStar;

    LineProducingSpecies(std::shared_ptr<Parameters> params) :
        parameters(params), quadrature(params), lambda(params){};

    void read(const Io& io, const Size l);
    void write(const Io& io, const Size l) const;

    void read_populations(const Io& io, const Size l, const string tag);
    void write_populations(const Io& io, const Size l, const string tag) const;

    inline Size index(const Size p, const Size i) const;

    inline Real get_emissivity(const Size p, const Size k) const;
    inline Real get_opacity(const Size p, const Size k) const;

    inline void check_for_convergence(const Real pop_prec);
    inline void check_for_convergence_trial(const Real pop_prec);

    inline void update_using_LTE(const Double2& abundance, const Vector<Real>& temperature);
    inline void set_LTE(
        const Double2& abundance, const Vector<Real>& temperature, const Size p, const Size k);

    inline void update_using_statistical_equilibrium(
        const Double2& abundance, const Vector<Real>& temperature);

    inline void update_using_statistical_equilibrium_sparse(
        const Double2& abundance, const Vector<Real>& temperature);

    inline void update_using_Ng_acceleration();
    inline void update_using_acceleration(const Size order);
    inline void update_using_acceleration_trial(const Size order);
    inline void correct_negative_populations(
        const Double2& abundance, const Vector<Real>& temperature);
};

#include "lineProducingSpecies.tpp"
