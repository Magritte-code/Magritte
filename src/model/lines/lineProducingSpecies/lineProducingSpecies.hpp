#pragma once

#include <Eigen/SparseLU>
using Eigen::SparseLU;
#include <Eigen/SparseCore>
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::Triplet;
using Eigen::COLAMDOrdering;


#include "io/io.hpp"
#include "tools/types.hpp"
#include "linedata/linedata.hpp"
#include "quadrature/quadrature.hpp"
#include "lambda/lambda.hpp"


class LineProducingSpecies
{
    public:
        Linedata   linedata;             ///< data for line producing species
        Quadrature quadrature;           ///< data for integral over line
        Lambda     lambda;               ///< Approximate Lambda Operator (ALO)

        Double2 Jlin;                    ///< actual mean intensity in the line
        Double2 Jeff;                    ///< effective mean intensity in the line (actual - ALO)
        Double2 Jdif;                    ///< effective mean intensity in the line (actual - ALO)

        Long3 nr_line;                   ///< frequency number corresponing to line (p,k,z)

        double relative_change_mean;     ///< mean    relative change
        double relative_change_max;      ///< maximum relative change
        double fraction_not_converged;   ///< fraction of levels that is not converged

        VectorXd population;             ///< level population (most recent)
        Double1  population_tot;         ///< total level population (sum over levels)

        VectorXd population_prev1;       ///< level populations 1 iteration  back
        VectorXd population_prev2;       ///< level populations 2 iterations back
        VectorXd population_prev3;       ///< level populations 3 iterations back

        SparseMatrix<double> RT;
        SparseMatrix<double> LambdaTest;
        SparseMatrix<double> LambdaStar;

        accel inline void set_npoints (const Size n);
        accel inline Size get_npoints () const;

        void read  (const Io& io, const Size l);
        void write (const Io& io, const Size l) const;

        void read_populations  (const Io &io, const Size l, const string tag);
        void write_populations (const Io &io, const Size l, const string tag) const;

        inline Size index (const Size p, const Size i) const;

        inline double get_emissivity (const Size p, const Size k) const;
        inline double get_opacity    (const Size p, const Size k) const;

        inline void check_for_convergence (
            const double pop_prec         );

        inline void update_using_LTE (
            const Double2 &abundance,
            const Double1 &temperature);

        inline void update_using_statistical_equilibrium (
            const Double2 &abundance,
            const Double1 &temperature                   );

        inline void update_using_Ng_acceleration ();

    private:
        Size nquads;
        Size npoints;
};


#include "lineProducingSpecies.tpp"