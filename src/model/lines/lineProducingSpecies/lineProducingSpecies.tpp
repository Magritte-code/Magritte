#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "tools/constants.hpp"
#include "tools/types.hpp"
#include "paracabs.hpp"

#include <map>



///  Indexer for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
///    @return corresponding index for p and i
//////////////////////////////////////////////
/// REQUIRED: Please do not change this unless you also have a look at solve_statistical_equilibrium for the 'compressed indices' and change them accordingly
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

    // populations.push_back (population);
}


inline void LineProducingSpecies :: check_for_convergence (const Real pop_prec, vector<Size> &points_in_grid)
{
    //Checking convergence for all points in the current grid
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

// DEPRECATED; no longer in use
// ///  update_using_acceleration: perform a Ng accelerated iteration step
// ///    for level populations. All variable names are based on lecture notes
// ///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
// ///////////////////////////////////////////////////////////////////////////
// void LineProducingSpecies :: update_using_acceleration (const Size order)
// {
//     MatrixXr RTR (order, order);
//
//     for (Size i = 0; i < order; i++)
//     {
//         for (Size j = 0; j < order; j++)
//         {
//             RTR(i,j) = residuals[i].dot(residuals[j]);
//         }
//     }
//
//     VectorXr ones  = VectorXr::Constant(order, 1.0);
//     VectorXr coef  = RTR.colPivHouseholderQr().solve(ones);
//              coef /= coef.sum();
//
//     residuals  .push_back(population-populations.back());
//     populations.push_back(population);
//
//     population = VectorXr::Zero(population.size());
//
//     for (Size i = 0; i < order; i++)
//     {
//         population += populations[order-1-i] * coef[order-1-i];
//     }
// }


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
    //, const Double2 levelpops_higher_level,
    //, Bool use_residual,
    //, const Double2 residualy_higher_level
{

    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    // residuals  .push_back(population-populations.back());
    // populations.push_back(population);

    std::cout<<"setting vector population"<<std::endl;
    VectorXr new_population = VectorXr::Zero (parameters.npoints()*linedata.nlev);

    const Size max_matrix_size=100*1024*1024;//hundred megabyte
    // const Size max_matrix_size=1*1024*1024;//one megabyte
    // It turns out that due to using triplets to construct the sparse matrix, the total extra memory needed for this calculation is about three times this number

    // const unsigned long long int total_non_zeros = parameters.npoints() * (      linedata.nlev
    //                                                + 6 * linedata.nrad
    //                                                + 4 * linedata.ncol_tot );

    // We can split our block diagonal matrix if there are no non-local contributions
    if (parameters.n_off_diag==0){
        Size n_points_per_block=max_matrix_size/((linedata.nlev + 6 * linedata.nrad + 4 * linedata.ncol_tot )*sizeof(Real));
        Size n_different_matrices=(points_in_grid.size()+(n_points_per_block-1))/n_points_per_block;//rounding up
        std::cout<<"Splitting big matrix into n parts: "<<n_different_matrices<<std::endl;

        //A trivial parallelisation would be possible, if it were not that we modify an internal state during the calculation
        for (Size idx=0;idx<n_different_matrices;idx++)
        {
            Size firstidx=idx*n_points_per_block;
            Size lastidx=std::min(static_cast<Size>(points_in_grid.size()-1),static_cast<Size>((idx+1)*n_points_per_block-1));
            vector<Size> current_points_in_block = std::vector<Size>(points_in_grid.begin() + firstidx, points_in_grid.begin()+lastidx+1);
            VectorXr resulting_y = solve_statistical_equilibrium(abundance,temperature,current_points_in_block);//does not need to be continguous

            //and finally putting those resulting y values at the right place in the y vector
            Size temp_idx=0;
            for (Size point_in_block:current_points_in_block)
            {
                new_population(Eigen::seq(index(point_in_block,0),index(point_in_block,linedata.nlev-1)))=resulting_y(Eigen::seq(index(temp_idx,0),index(temp_idx,linedata.nlev-1)));
                temp_idx++;
            }
        }
    }
    else
    {//unfortunately, we cannot split this matrix into smaller pieces without making any errors
      new_population = solve_statistical_equilibrium(abundance,temperature,points_in_grid);
    }

    cout << "Succesfully solved for the level populations!"       << endl;

    population=new_population;

}


///  Refactor of the matrix solving for allowing to use less point per block matrix (decreases memory usage)
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
///    @param[in] points_in_grid: the points in the current grid
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// IMPORTANT: Only call this with a part of the points currently in grid if parameters.n_off_diag==0
/// OTHERWISE this will result in bogus (starting with the approximated lambda operator part)
inline VectorXr LineProducingSpecies::solve_statistical_equilibrium(const Double2 &abundance,
        const Vector<Real> &temperature, vector<Size> &points_to_use)
{

    if ((parameters.n_off_diag!=0)&&(parameters.npoints()!=points_to_use.size()))//if you do not set n_off_diag to 0, then we cannot split this almost 'block'matrix into blocks
    {// if parameters.n_off_diag!=0, then the points are assumed to be ordered from low to high
        std::cout<<"Attempted to split non-block matrix into blocks"<<std::endl;
        throw std::runtime_error("Attempted to split non-block matrix into blocks");
    }


    VectorXr y = VectorXr::Zero (points_to_use.size()*linedata.nlev);


    const Size non_zeros = points_to_use.size() * (      linedata.nlev
                                                     + 6 * linedata.nrad
                                                     + 4 * linedata.ncol_tot );


    vector<Triplet<Real, Size>> triplets;

    triplets   .reserve (non_zeros);

    Size nbpoints=points_to_use.size();

    std::map<Size,Size> indexmap;//for simplicity, maps the indices from points_to_use to consecutive indices
    Size temp_idx=0;
    for(Size point: points_to_use)
    {
        indexmap.insert(std::pair<Size,Size>(point,temp_idx));
        temp_idx++;
    }


    for (Size p : points_to_use)
    {
        // Radiative transitions

        for (Size k = 0; k < linedata.nrad; k++)
          {
              const Real v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
              const Real v_JI =                 linedata.Ba[k] * Jeff[p][k];

              const Real t_IJ = linedata.Bs[k] * Jdif[p][k];
              const Real t_JI = linedata.Ba[k] * Jdif[p][k];

              // Note: we define our transition matrix as the transpose of R in the paper.
              // const Size I = index (p, linedata.irad[k]);
              // const Size J = index (p, linedata.jrad[k]);

              //compressed index notation (such that there are no zero rows when points_to_use!=all grid points)
              const Size I = index (indexmap.at(p), linedata.irad[k]);
              const Size J = index (indexmap.at(p), linedata.jrad[k]);


              if (linedata.jrad[k] != linedata.nlev-1)
              {

                  triplets   .push_back (Triplet<Real, Size> (J, I, +v_IJ));
                  triplets   .push_back (Triplet<Real, Size> (J, J, -v_JI));

                  // triplets_LS.push_back (Triplet<Real, Size> (J, I, +t_IJ));
                  // triplets_LS.push_back (Triplet<Real, Size> (J, J, -t_JI));
              }

              if (linedata.irad[k] != linedata.nlev-1)
              {

                  triplets   .push_back (Triplet<Real, Size> (I, J, +v_JI));
                  triplets   .push_back (Triplet<Real, Size> (I, I, -v_IJ));

                  // triplets_LS.push_back (Triplet<Real, Size> (I, J, +t_JI));
                  // triplets_LS.push_back (Triplet<Real, Size> (I, I, -t_IJ));
              }
          }

          // Approximated Lambda operator

          for (Size k = 0; k < linedata.nrad; k++)
          {
              for (Size m = 0; m < lambda.get_size(p,k); m++)
              {
                  const Size   nr =  lambda.get_nr(p, k, m);
                  const Real v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                  // Note: we define our transition matrix as the transpose of R in the paper.
                  // Compressed index notation (such that there are no zero rows when points_to_use!=all grid points)
                  const Size I = index (indexmap.at(nr), linedata.irad[k]);
                  const Size J = index (indexmap.at(p),  linedata.jrad[k]);

                  if (linedata.jrad[k] != linedata.nlev-1)
                  {
                      triplets   .push_back (Triplet<Real, Size> (J, I, +v_IJ));
                      // triplets_LT.push_back (Triplet<Real, Size> (J, I, +v_IJ));
                  }

                  if (linedata.irad[k] != linedata.nlev-1)
                  {
                      triplets   .push_back (Triplet<Real, Size> (I, I, -v_IJ));
                      // triplets_LT.push_back (Triplet<Real, Size> (I, I, -v_IJ));
                  }
              }
          }

          // Collisional transitions

          for (CollisionPartner &colpar : linedata.colpar)
          {
              Real abn = abundance[p][colpar.num_col_partner];
              Real tmp = temperature[p];

              colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
              colpar.interpolate_collision_coefficients (tmp);//CURRENTLY NOT THREAD SAFE
              //ADJUSTS INTERNAL STATE

              for (Size k = 0; k < colpar.ncol; k++)
              {
                  const Real v_IJ = colpar.Cd_intpld[k] * abn;//ONLY HERE IS THIS INTERNAL STATE USED
                  const Real v_JI = colpar.Ce_intpld[k] * abn;

                  // Note: we define our transition matrix as the transpose of R in the paper.
                  // const Size I = index (p, colpar.icol[k]);
                  // const Size J = index (p, colpar.jcol[k]);

                  // Compressed index notation (such that there are no zero rows when points_to_use!=all grid points)
                  const Size I = index (indexmap.at(p), colpar.icol[k]);
                  const Size J = index (indexmap.at(p), colpar.jcol[k]);


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
              // const Size I = index (p, linedata.nlev-1);
              // const Size J = index (p, i);

              //compressed index notation (such that there are no zero rows when points_to_use!=all grid points)
              const Size I = index (indexmap.at(p), linedata.nlev-1);
              const Size J = index (indexmap.at(p), i);

              triplets.push_back (Triplet<Real, Size> (I, J, 1.0));
          }

          // y[index (p, linedata.nlev-1)] = population_tot[p];
          // Again compressed index notation
          y[index (indexmap.at(p), linedata.nlev-1)] = population_tot[p];

      } // for all cells

      // std::cout<<"Setting RT"<<std::endl;

      Eigen::SparseMatrix<Real> RTmat;

      RTmat.resize (points_to_use.size()*linedata.nlev,
                    points_to_use.size()*linedata.nlev );

      RTmat.setFromTriplets (triplets.begin(), triplets.end());

      // std::cout << "The matrix RT is of size "
      //     << RTmat.rows() << "x" << RTmat.cols() << std::endl;
      //
      // std::cout << "The vector y is of size " << y.size() << std::endl;

      SparseLU <SparseMatrix<Real>, COLAMDOrdering<int>> solver;

      // cout << "Analyzing system of rate equations..."      << endl;

      solver.analyzePattern (RTmat);

      // cout << "Factorizing system of rate equations..."    << endl;

      solver.factorize (RTmat);

      if (solver.info() != Eigen::Success)
      {
          cout << "Factorization failed with error message:" << endl;
          cout << solver.lastErrorMessage()                  << endl;
          // cout << endl << RT << endl;

          throw std::runtime_error ("Eigen solver ERROR.");
      }

      // cout << "Solving rate equations for the level populations..." << endl;

      VectorXr temp_population = solver.solve (y);

      if (solver.info() != Eigen::Success)
      {
          cout << "Solving failed with error:" << endl;
          cout << solver.lastErrorMessage()    << endl;
          assert (false);
      }

      // cout << "Succesfully solved for the level populations!"       << endl;

      return temp_population;

}

///  Sets all the level populations of this species to the new_population
///    @param[in] new_population: The new level populations of this species
inline void LineProducingSpecies::set_all_level_pops(VectorXr new_population)
{
  population=new_population;
}

///  Returns all level populations of this species
///    @return The current level populations of this species
inline VectorXr LineProducingSpecies::get_all_level_pops()
{
  return population;
}
