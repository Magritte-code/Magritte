#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "tools/constants.hpp"
#include "tools/types.hpp"
#include "paracabs.hpp"


///  Indexer for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
///    @return corresponding index for p and i
//////////////////////////////////////////////
inline Size LineProducingSpecies :: index (const Size p, const Size i) const
{
    return i + p*linedata.nlev;
}
/// Getter for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
inline Real LineProducingSpecies :: get_level_pop  (const Size p, const Size i) const
{
    return population(index(p,i));
}
/// Setter for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
inline void LineProducingSpecies :: set_level_pop  (const Size p, const Size i, const Real value)
{
    population(index(p,i))=value;
}


///  Getter for the line emissivity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line emissivity for cell p and transition k
//////////////////////////////////////////////////////////
inline Real LineProducingSpecies :: get_emissivity (const Size p, const Size k) const
{
  const Size i = index (p, linedata.irad[k]);

  return HH_OVER_FOUR_PI * linedata.A[k] * population(i);
}


///  Getter for the line opacity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line opacity for cell p and transition k
///////////////////////////////////////////////////////
inline Real LineProducingSpecies :: get_opacity (const Size p, const Size k) const
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
    const Double2      &abundance,
    const Vector<Real> &temperature )
{
    threaded_for (p, parameters.npoints(),
    {
        population_tot[p] = abundance[p][linedata.num];

        Real partition_function = 0.0;

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const Size ind = index (p, i);

            population(ind) = linedata.weight[i]
                              * exp (-linedata.energy[i] / (KB*temperature[p]));

            partition_function += population (ind);
        }

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const Size ind = index (p, i);

            population(ind) *= population_tot[p] / partition_function;
        }
    })

    populations.push_back (population);
}


inline void LineProducingSpecies :: check_for_convergence (const Real pop_prec, vector<Size> &points_in_grid)
{
    //this also needs some changes for coarser grids; because interpolation happens at later time, there can be inconsistencies
    const Size nbpoints=points_in_grid.size();
    const Real weight = 1.0 / (nbpoints * linedata.nlev);
    // const Real weight = 1.0 / (parameters.npoints() * linedata.nlev);

    Real fnc = 0.0;
    Real rcm = 0.0;

    relative_change_max = 0.0;

//    for (long p = 0; p < ncells; p++)
#   pragma omp parallel for reduction (+: fnc, rcm)
    for (Size idx = 0; idx < nbpoints; idx++)
    {
        const Size p=points_in_grid[idx];
        const double min_pop = 1.0E-10 * population_tot[p];

        for (Size i = 0; i < linedata.nlev; i++)
        {
            const Size ind = index (p, i);

            if (population(ind) > min_pop)
            {
                Real relative_change = 2.0;

                //note: with multigrid, the population at points not in the grid does not get updated quickly enough, making the relative change always 0
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
    VectorXr Wt (parameters.npoints()*linedata.nlev);

    VectorXr Q1 = population - 2.0*population_prev1 + population_prev2;
    VectorXr Q2 = population -     population_prev1 - population_prev2 + population_prev3;
    VectorXr Q3 = population -     population_prev1;

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

    const Real A1 = Q1.dot(Q1);
    const Real A2 = Q1.dot(Q2);
    const Real B2 = Q2.dot(Q2);
    const Real C1 = Q1.dot(Q3);
    const Real C2 = Q2.dot(Q3);

    const Real B1 = A2;

    const Real denominator = A1*B2 - A2*B1;

    if (denominator != 0.0)
    {
        const VectorXr pop_tmp = population;

        const Real a = (C1*B2 - C2*B1) / denominator;
        const Real b = (C2*A1 - C1*A2) / denominator;

        population = (1.0 - a - b)*population
                               + a*population_prev1
                               + b*population_prev2;

        population_prev3 = population_prev2;
        population_prev2 = population_prev1;
        population_prev1 = pop_tmp;
    }
}


///  update_using_acceleration: perform a Ng accelerated iteration step
///    for level populations. All variable names are based on lecture notes
///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
///////////////////////////////////////////////////////////////////////////
void LineProducingSpecies :: update_using_acceleration (const Size order)
{
    MatrixXr RTR (order, order);

    for (Size i = 0; i < order; i++)
    {
        for (Size j = 0; j < order; j++)
        {
            RTR(i,j) = residuals[i].dot(residuals[j]);
        }
    }

    VectorXr ones  = VectorXr::Constant(order, 1.0);
    VectorXr coef  = RTR.colPivHouseholderQr().solve(ones);
             coef /= coef.sum();

    residuals  .push_back(population-populations.back());
    populations.push_back(population);

    population = VectorXr::Zero(population.size());

    for (Size i = 0; i < order; i++)
    {
        population += populations[order-1-i] * coef[order-1-i];
    }
}


///  update_using_statistical_equilibrium: computes level populations by solving
///  the statistical equilibrium equation taking into account the radiation field
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
///    @param[in] points_in_grid: the points in the current grid
/////////////////////////////////////////////////////////////////////////////////
/// Note: for points currently not in the grid, the solution may be bogus
/// FIXME: instead return the previous values for the level pops
inline void LineProducingSpecies :: update_using_statistical_equilibrium (
    const Double2      &abundance,
    const Vector<Real> &temperature,
    vector<Size> &points_in_grid)
{

    const Size non_zeros = parameters.npoints() * (      linedata.nlev
                                                   + 6 * linedata.nrad
                                                   + 4 * linedata.ncol_tot );

    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    residuals  .push_back(population-populations.back());
    populations.push_back(population);

//    SparseMatrix<double> RT (ncells*linedata.nlev, ncells*linedata.nlev);

    VectorXr y = VectorXr::Zero (parameters.npoints()*linedata.nlev);

    vector<Triplet<Real, Size>> triplets;
    vector<Triplet<Real, Size>> triplets_LT;
    vector<Triplet<Real, Size>> triplets_LS;

    triplets   .reserve (non_zeros);
    triplets_LT.reserve (non_zeros);
    triplets_LS.reserve (non_zeros);

    //todo: use current grid ...
    // vector<Size> points_in_grid=model.geometry.points.multiscale.get_current_points_in_grid();
    Size nbpoints=points_in_grid.size();

    for (Size idx = 0; idx < nbpoints; idx++) // !!! no OMP because push_back is not thread safe !!!
    {
        const Size p=points_in_grid[idx];
        std::cout<<"current point: "<<p<<std::endl;
        // Radiative transitions

        for (Size k = 0; k < linedata.nrad; k++)
        {
            const Real v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
            const Real v_JI =                 linedata.Ba[k] * Jeff[p][k];

            const Real t_IJ = linedata.Bs[k] * Jdif[p][k];
            const Real t_JI = linedata.Ba[k] * Jdif[p][k];

            // Note: we define our transition matrix as the transpose of R in the paper.
            const Size I = index (p, linedata.irad[k]);
            const Size J = index (p, linedata.jrad[k]);

            std::cout<<"I: "<<I<<std::endl;
            std::cout<<"J: "<<J<<std::endl;

            if (linedata.jrad[k] != linedata.nlev-1)
            {
                triplets   .push_back (Triplet<Real, Size> (J, I, +v_IJ));
                triplets   .push_back (Triplet<Real, Size> (J, J, -v_JI));

                triplets_LS.push_back (Triplet<Real, Size> (J, I, +t_IJ));
                triplets_LS.push_back (Triplet<Real, Size> (J, J, -t_JI));
            }

            if (linedata.irad[k] != linedata.nlev-1)
            {
                triplets   .push_back (Triplet<Real, Size> (I, J, +v_JI));
                triplets   .push_back (Triplet<Real, Size> (I, I, -v_IJ));

                triplets_LS.push_back (Triplet<Real, Size> (I, J, +t_JI));
                triplets_LS.push_back (Triplet<Real, Size> (I, I, -t_IJ));
            }
        }

        std::cout<<"flag 2"<<std::endl;

        // Approximated Lambda operator

        for (Size k = 0; k < linedata.nrad; k++)
        {
            for (Size m = 0; m < lambda.get_size(p,k); m++)
            {
                const Size   nr =  lambda.get_nr(p, k, m);
                const Real v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                // Note: we define our transition matrix as the transpose of R in the paper.
                const Size I = index (nr, linedata.irad[k]);
                const Size J = index (p,  linedata.jrad[k]);

                std::cout<<"I: "<<I<<std::endl;
                std::cout<<"J: "<<J<<std::endl;

                if (linedata.jrad[k] != linedata.nlev-1)
                {
                    triplets   .push_back (Triplet<Real, Size> (J, I, +v_IJ));
                    triplets_LT.push_back (Triplet<Real, Size> (J, I, +v_IJ));
                }

                if (linedata.irad[k] != linedata.nlev-1)
                {
                    triplets   .push_back (Triplet<Real, Size> (I, I, -v_IJ));
                    triplets_LT.push_back (Triplet<Real, Size> (I, I, -v_IJ));
                }
            }
        }

        std::cout<<"flag 3"<<std::endl;

        // Collisional transitions

        for (CollisionPartner &colpar : linedata.colpar)
        {
            Real abn = abundance[p][colpar.num_col_partner];
            Real tmp = temperature[p];

            colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
            colpar.interpolate_collision_coefficients (tmp);


            for (Size k = 0; k < colpar.ncol; k++)
            {
                const Real v_IJ = colpar.Cd_intpld[k] * abn;
                const Real v_JI = colpar.Ce_intpld[k] * abn;


                // Note: we define our transition matrix as the transpose of R in the paper.
                const Size I = index (p, colpar.icol[k]);
                const Size J = index (p, colpar.jcol[k]);

                std::cout<<"I: "<<I<<std::endl;
                std::cout<<"J: "<<J<<std::endl;

                if (colpar.jcol[k] != linedata.nlev-1)
                {
                    triplets.push_back (Triplet<Real, Size> (J, I, +v_IJ));
                    triplets.push_back (Triplet<Real, Size> (J, J, -v_JI));
                }

                if (colpar.icol[k] != linedata.nlev-1)
                {
                    triplets.push_back (Triplet<Real, Size> (I, J, +v_JI));
                    triplets.push_back (Triplet<Real, Size> (I, I, -v_IJ));
                }
            }
        }


        for (Size i = 0; i < linedata.nlev; i++)
        {
            const Size I = index (p, linedata.nlev-1);
            const Size J = index (p, i);

            std::cout<<"I: "<<I<<std::endl;
            std::cout<<"J: "<<J<<std::endl;

            triplets.push_back (Triplet<Real, Size> (I, J, 1.0));
        }

        //y.insert (index (p, linedata.nlev-1)) = population_tot[p];
        y[index (p, linedata.nlev-1)] = population_tot[p];
        //y.insert (index (p, linedata.nlev-1)) = 1.0;//population_tot[p];


    } // for all cells


    // adding 1 to diagonal for line transitions currently not in the grid
    // currently, they should not be affected in any way (only interpolation may change them)

    // well, points_in_grid should be ordered, but extra precautions can never hurt
    std::sort(points_in_grid.begin(),points_in_grid.end());
    Size otheridx=0;
    for (Size p=0;p<parameters.npoints();p++)
    {
      if (p==points_in_grid[otheridx])
      {otheridx++;}
      else
      {//add 1 to all diagonal line transitions of this point
        for (Size l = 0; l < linedata.nlev; l++)
        {
          const Size I = index (p, l);
          triplets.push_back (Triplet<Real, Size> (I, I, 1.0));
          // RT.insert(I,I)=1.0;
          std::cout<<"added to index: "<<I<<std::endl;
        }
      }
    }

    RT        .setFromTriplets (triplets   .begin(), triplets   .end());
    LambdaStar.setFromTriplets (triplets_LS.begin(), triplets_LS.end());
    LambdaTest.setFromTriplets (triplets_LT.begin(), triplets_LT.end());

    std::cout<<"RT rows: "<<RT.rows()<<std::endl;
    std::cout<<"RT cols: "<<RT.cols()<<std::endl;





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


    SparseLU <SparseMatrix<Real>, COLAMDOrdering<int>> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;


    cout << "Analyzing system of rate equations..."      << endl;

    solver.analyzePattern (RT);

    cout << "Factorizing system of rate equations..."    << endl;

    solver.factorize (RT);

    if (solver.info() != Eigen::Success)
    {
        cout << "Factorization failed with error message:" << endl;
        cout << solver.lastErrorMessage()                  << endl;

        // cout << endl << RT << endl;

        throw std::runtime_error ("Eigen solver ERROR.");
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

    // VectorXr reduced_population = solver.solve (y);
    //
    // for (Size p: points_in_grid)
    // {
    //   population(p) = reduced_population(p);
    // }

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
