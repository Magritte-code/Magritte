#include <math.h>
#include <Eigen/Core>

#include "tools/constants.hpp"
#include "tools/types.hpp"


accel inline void LineProducingSpecies :: set_npoints (const Size n)
{
    npoints = n;
}


accel inline Size LineProducingSpecies :: get_npoints () const
{
    return npoints;
}


///  Indexer for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
///    @return corresponding index for p and i
//////////////////////////////////////////////
inline Size LineProducingSpecies :: index (const Size p, const Size i) const
{
    return i + p*linedata.nlev;
}


///  Getter for the line emissivity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line emissivity for cell p and transition k
//////////////////////////////////////////////////////////
inline double LineProducingSpecies :: get_emissivity (const Size p, const Size k) const
{
  const Size i = index (p, linedata.irad[k]);

  return HH_OVER_FOUR_PI * linedata.A[k] * population(i);
}


///  Getter for the line opacity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line opacity for cell p and transition k
///////////////////////////////////////////////////////
inline double LineProducingSpecies :: get_opacity (const Size p, const Size k) const
{
  const Size i = index (p, linedata.irad[k]);
  const Size j = index (p, linedata.jrad[k]);

  return HH_OVER_FOUR_PI * (  population(j) * linedata.Ba[k]
                            - population(i) * linedata.Bs[k] );
}


///  set_LTE_level_populations
///    @param[in] abundance_lspec: abundance of line species
///    @param[in] temperature: local gas temperature
///    @param[in] p: number of cell
///    @param[in] l: number of line producing species
///////////////////////////////////////////////////////////
inline void LineProducingSpecies :: update_using_LTE (
    const Double2 &abundance,
    const Double1 &temperature )
{
    threaded_for (p, ncells)
    {
        population_tot[p] = abundance[p][linedata.num];

        double partition_function = 0.0;

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const long ind = index (p, i);

            population(ind) = linedata.weight[i]
                              * exp (-linedata.energy[i] / (KB*temperature[p]));

            partition_function += population (ind);
        }

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const long ind = index (p, i);

            population(ind) *= population_tot[p] / partition_function;
        }
    }
}




inline void LineProducingSpecies :: check_for_convergence (const double pop_prec)
{
    const double weight = 1.0 / (ncells * linedata.nlev);

    double fnc = 0.0;
    double rcm = 0.0;

    relative_change_max = 0.0;


//    for (long p = 0; p < ncells; p++)
#   pragma omp parallel for reduction (+: fnc, rcm)
    for (Size p = 0; p < ncells; p++)
    {
        const double min_pop = 1.0E-10 * population_tot[p];

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const long ind = index (p, i);

            if (population(ind) > min_pop)
            {
                double relative_change = 2.0;

                relative_change *= fabs (population (ind) - population_prev1 (ind));
                relative_change /=      (population (ind) + population_prev1 (ind));

                if (relative_change > pop_prec)
                {
                    fnc += weight;
                }

                rcm += (weight * relative_change);

                // NOT THREAD SAFE !!!
                //if (relative_change > relative_change_max)
                //{
                //  relative_change_max = relative_change;
                //}
            }
        }
    }

    fraction_not_converged = fnc;
    relative_change_mean   = rcm;
}




///  update_using_Ng_acceleration: perform a Ng accelerated iteration step
///    for level populations. All variable names are based on lecture notes
///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
///////////////////////////////////////////////////////////////////////////
void LineProducingSpecies :: update_using_Ng_acceleration ()
{
    VectorXd Wt (ncells*linedata.nlev);

    VectorXd Q1 = population - 2.0*population_prev1 + population_prev2;
    VectorXd Q2 = population -     population_prev1 - population_prev2 + population_prev3;
    VectorXd Q3 = population -     population_prev1;

    //OMP_PARALLEL_FOR (ind, ncells*linedata.nlev)
    //{
    //  if (population (ind) > 0.0)
    //  {
    //    Wt (ind) = Jlin[p][k];
    //  }

    //  else
    //  {
    //    Wt (ind) = 1.0;
    //  }
    //}

    //const double A1 = Q1.dot (Wt.asDiagonal()*Q1);
    //const double A2 = Q1.dot (Wt.asDiagonal()*Q2);
    //const double B2 = Q2.dot (Wt.asDiagonal()*Q2);
    //const double C1 = Q1.dot (Wt.asDiagonal()*Q3);
    //const double C2 = Q2.dot (Wt.asDiagonal()*Q3);

    const double A1 = Q1.dot(Q1);
    const double A2 = Q1.dot(Q2);
    const double B2 = Q2.dot(Q2);
    const double C1 = Q1.dot(Q3);
    const double C2 = Q2.dot(Q3);

    const double B1 = A2;

    const double denominator = A1*B2 - A2*B1;


    if (denominator != 0.0)
    {
        const VectorXd pop_tmp = population;

        const double a = (C1*B2 - C2*B1) / denominator;
        const double b = (C2*A1 - C1*A2) / denominator;

        population = (1.0 - a - b)*population
                               + a*population_prev1
                               + b*population_prev2;

        population_prev3 = population_prev2;
        population_prev2 = population_prev1;
        population_prev1 = pop_tmp;
    }

}




///  update_using_statistical_equilibrium: computes level populations by solving
///  the statistical equilibrium equation taking into account the radiation field
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
/////////////////////////////////////////////////////////////////////////////////

inline void LineProducingSpecies ::
    update_using_statistical_equilibrium (
        const Double2 &abundance,
        const Double1 &temperature       )
{

    const size_t non_zeros = ncells * (      linedata.nlev
                                       + 6 * linedata.nrad
                                       + 4 * linedata.ncol_tot );


    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;


//    SparseMatrix<double> RT (ncells*linedata.nlev, ncells*linedata.nlev);

    VectorXd y = VectorXd::Zero (ncells*linedata.nlev);

    vector<Triplet<double, int>> triplets;
    vector<Triplet<double, int>> triplets_LT;
    vector<Triplet<double, int>> triplets_LS;

    triplets   .reserve (non_zeros);
    triplets_LT.reserve (non_zeros);
    triplets_LS.reserve (non_zeros);

    for (size_t p = 0; p < ncells; p++) // !!! no OMP because push_back is not thread safe !!!
    {

        // Radiative transitions

        for (size_t k = 0; k < linedata.nrad; k++)
        {
            const double v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
            const double v_JI =                 linedata.Ba[k] * Jeff[p][k];

            const double t_IJ = linedata.Bs[k] * Jdif[p][k];
            const double t_JI = linedata.Ba[k] * Jdif[p][k];


            // Note: we define our transition matrix as the transpose of R in the paper.
            const size_t I = index (p, linedata.irad[k]);
            const size_t J = index (p, linedata.jrad[k]);

            if (linedata.jrad[k] != linedata.nlev-1)
            {
                triplets   .push_back (Triplet<double, int> (J, I, +v_IJ));
                triplets   .push_back (Triplet<double, int> (J, J, -v_JI));

                triplets_LS.push_back (Triplet<double, int> (J, I, +t_IJ));
                triplets_LS.push_back (Triplet<double, int> (J, J, -t_JI));
            }

            if (linedata.irad[k] != linedata.nlev-1)
            {
                triplets   .push_back (Triplet<double, int> (I, J, +v_JI));
                triplets   .push_back (Triplet<double, int> (I, I, -v_IJ));

                triplets_LS.push_back (Triplet<double, int> (I, J, +t_JI));
                triplets_LS.push_back (Triplet<double, int> (I, I, -t_IJ));
            }
        }


        // Approximated Lambda operator

        for (size_t k = 0; k < linedata.nrad; k++)
        {
            for (size_t m = 0; m < lambda.get_size(p,k); m++)
            {
                const size_t   nr =  lambda.get_nr(p, k, m);
                const double v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                // Note: we define our transition matrix as the transpose of R in the paper.
                const size_t I = index (nr, linedata.irad[k]);
                const size_t J = index (p,  linedata.jrad[k]);

                if (linedata.jrad[k] != linedata.nlev-1)
                {
                    triplets   .push_back (Triplet<double, int> (J, I, +v_IJ));
                    triplets_LT.push_back (Triplet<double, int> (J, I, +v_IJ));
                }

                if (linedata.irad[k] != linedata.nlev-1)
                {
                    triplets   .push_back (Triplet<double, int> (I, I, -v_IJ));
                    triplets_LT.push_back (Triplet<double, int> (I, I, -v_IJ));
                }
            }
        }



        // Collisional transitions

        for (CollisionPartner &colpar : linedata.colpar)
        {
            double abn = abundance[p][colpar.num_col_partner];
            double tmp = temperature[p];

            colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
            colpar.interpolate_collision_coefficients (tmp);


            for (size_t k = 0; k < colpar.ncol; k++)
            {
                const double v_IJ = colpar.Cd_intpld[k] * abn;
                const double v_JI = colpar.Ce_intpld[k] * abn;


                // Note: we define our transition matrix as the transpose of R in the paper.
                const size_t I = index (p, colpar.icol[k]);
                const size_t J = index (p, colpar.jcol[k]);

                if (colpar.jcol[k] != linedata.nlev-1)
                {
                    triplets.push_back (Triplet<double, int> (J, I, +v_IJ));
                    triplets.push_back (Triplet<double, int> (J, J, -v_JI));
                }

                if (colpar.icol[k] != linedata.nlev-1)
                {
                    triplets.push_back (Triplet<double, int> (I, J, +v_JI));
                    triplets.push_back (Triplet<double, int> (I, I, -v_IJ));
                }
            }
        }


        for (size_t i = 0; i < linedata.nlev; i++)
        {
            const long I = index (p, linedata.nlev-1);
            const long J = index (p, i);

            triplets.push_back (Triplet<double, int> (I, J, 1.0));
        }

        //y.insert (index (p, linedata.nlev-1)) = population_tot[p];
        y[index (p, linedata.nlev-1)] = population_tot[p];
        //y.insert (index (p, linedata.nlev-1)) = 1.0;//population_tot[p];

    } // for all cells


    RT        .setFromTriplets (triplets   .begin(), triplets   .end());
    LambdaStar.setFromTriplets (triplets_LS.begin(), triplets_LS.end());
    LambdaTest.setFromTriplets (triplets_LT.begin(), triplets_LT.end());


    //cout << "Compressing RT" << endl;

    //RT.makeCompressed ();


    //Eigen::BiCGSTAB <SparseMatrix<double>> solver;

    //cout << "Try compute" << endl;

    //solver.compute (RT);

    //if (solver.info() != Eigen::Success)
    //{
    //  cout << "Decomposition failed" << endl;
    //  //assert(false);
    //}


    //for (int tel=0; tel<5; tel++)
    //{
    //  //Eigen::Gues x0 = population;

    //  population = solver.solveWithGuess (y, population);
    //  std::cout << "#iterations:     " << solver.iterations() << std::endl;
    //  std::cout << "estimated error: " << solver.error()      << std::endl;
    //}

    //assert (false);


    SparseLU <SparseMatrix<double>, COLAMDOrdering<int>> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;


    cout << "Analyzing system of rate equations..."      << endl;

    solver.analyzePattern (RT);

    cout << "Factorizing system of rate equations..."    << endl;

    solver.factorize (RT);

    if (solver.info() != Eigen::Success)
    {
        cout << "Factorization failed with error message:" << endl;
        cout << solver.lastErrorMessage()                  << endl;

        cout << endl << RT << endl;

        assert (false);
    }


    //cout << "Try compute" << endl;

    //solver.compute (RT);

    //if (solver.info() != Eigen::Success)
    //{
    //  cout << "Decomposition failed" << endl;
    //  //assert(false);
    //}


    cout << "Solving rate equations for the level populations..." << endl;

    population = solver.solve (y);

    if (solver.info() != Eigen::Success)
    {
        cout << "Solving failed with error:" << endl;
        cout << solver.lastErrorMessage()    << endl;
        assert (false);
    }

    cout << "Succesfully solved for the level populations!"       << endl;


    //OMP_PARALLEL_FOR (p, ncells)
    //{
    //
    //  for (long i = 0; i < linedata.nlev; i++)
    //  {
    //    const long I = index (p, i);

    //    population[I] = population_prev1[I];

    //    //if (population[I] < 1.0E-50)
    //    //{
    //    //  population[I] = 1.0E-50;
    //    //}
    //  }
    //}


}
