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
    threaded_for (p, parameters->npoints(),
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


inline void LineProducingSpecies :: check_for_convergence (const Real pop_prec)
{
    const Real weight = 1.0 / (parameters->npoints() * linedata.nlev);

    Real fnc = 0.0;
    Real rcm = 0.0;

#   pragma omp parallel for reduction (+: fnc, rcm)
    for (Size p = 0; p < parameters->npoints(); p++)
    {
        for (Size i = 0; i < linedata.nlev; i++)
        {
            const Size ind = index (p, i);

            if (population(ind) > parameters->min_rel_pop_for_convergence * population_tot[p])
            {
                Real relative_change = 2.0;

                relative_change *= fabs (population (ind) - population_prev1 (ind));
                relative_change /=      (population (ind) + population_prev1 (ind));

                if (relative_change > pop_prec)
                {
                    fnc += weight;
                }

                rcm += (weight * relative_change);
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
    VectorXr Wt (parameters->npoints()*linedata.nlev);

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
/////////////////////////////////////////////////////////////////////////////////
inline void LineProducingSpecies :: update_using_statistical_equilibrium (
    const Double2      &abundance,
    const Vector<Real> &temperature )
{
    RT        .resize (parameters->npoints()*linedata.nlev, parameters->npoints()*linedata.nlev);
    LambdaStar.resize (parameters->npoints()*linedata.nlev, parameters->npoints()*linedata.nlev);
    LambdaTest.resize (parameters->npoints()*linedata.nlev, parameters->npoints()*linedata.nlev);

    const Size non_zeros = parameters->npoints() * (      linedata.nlev
                                                    + 6 * linedata.nrad
                                                    + 4 * linedata.ncol_tot );
    // Store previous iterations
    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    residuals  .push_back(population-populations.back());
    populations.push_back(population);

//    SparseMatrix<double> RT (ncells*linedata.nlev, ncells*linedata.nlev);

    VectorXr y = VectorXr::Zero (parameters->npoints()*linedata.nlev);

    vector<Triplet<Real, Size>> triplets;
//    vector<Triplet<Real, Size>> triplets_LT;
//    vector<Triplet<Real, Size>> triplets_LS;

    triplets   .reserve (non_zeros);
//    triplets_LT.reserve (non_zeros);
//    triplets_LS.reserve (non_zeros);

    for (Size p = 0; p < parameters->npoints(); p++) // !!! no OMP because push_back is not thread safe !!!
    {
        // Radiative transitions

        for (Size k = 0; k < linedata.nrad; k++)
        {
            const Real v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
            const Real v_JI =                 linedata.Ba[k] * Jeff[p][k];

            // const Real t_IJ = linedata.Bs[k] * Jdif[p][k];
            // const Real t_JI = linedata.Ba[k] * Jdif[p][k];

            // Note: we define our transition matrix as the transpose of R in the paper.
            const Size I = index (p, linedata.irad[k]);
            const Size J = index (p, linedata.jrad[k]);

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
                const Size I = index (nr, linedata.irad[k]);
                const Size J = index (p,  linedata.jrad[k]);

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
            colpar.interpolate_collision_coefficients (tmp);

            for (Size k = 0; k < colpar.ncol; k++)
            {
                const Real v_IJ = colpar.Cd_intpld[k] * abn;
                const Real v_JI = colpar.Ce_intpld[k] * abn;

                // Note: we define our transition matrix as the transpose of R in the paper.
                const Size I = index (p, colpar.icol[k]);
                const Size J = index (p, colpar.jcol[k]);

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

            triplets.push_back (Triplet<Real, Size> (I, J, 1.0));
        }

        y[index (p, linedata.nlev-1)] = population_tot[p];

    } // for all cells


    RT        .setFromTriplets (triplets   .begin(), triplets   .end());
    // LambdaStar.setFromTriplets (triplets_LS.begin(), triplets_LS.end());
    // LambdaTest.setFromTriplets (triplets_LT.begin(), triplets_LT.end());


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

        throw std::runtime_error ("Eigen solver ERROR.");
    }

    cout << "Solving rate equations for the level populations..." << endl;

    population = solver.solve (y);

    if (solver.info() != Eigen::Success)
    {
        cout << "Solving failed with error:" << endl;
        cout << solver.lastErrorMessage()    << endl;
        assert (false);
    }

    cout << "Succesfully solved for the level populations!"       << endl;

}


///  update_using_statistical_equilibrium: computes level populations by solving
///  the statistical equilibrium equation taking into account the radiation field
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
/////////////////////////////////////////////////////////////////////////////////
inline void LineProducingSpecies :: update_using_statistical_equilibrium_sparse (
    const Double2      &abundance,
    const Vector<Real> &temperature )
{

    if (parameters->n_off_diag != 0)
    {
        throw std::runtime_error ("parameters->n_off_diag != 0 so cannot use sparse update.");
    }

    // Store previous iterations
    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    // residuals  .push_back(population-populations.back());
    // populations.push_back(population);

    pc::multi_threading::ThreadPrivate<MatrixXr> StatEq_;
    pc::multi_threading::ThreadPrivate<VectorXr>      y_;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        StatEq_(i).resize(linedata.nlev, linedata.nlev);

        y_(i) = VectorXr::Zero(linedata.nlev);
    }


    threaded_for (p, parameters->npoints(),
    {
        MatrixXr& StatEq = StatEq_();
        VectorXr&      y =      y_();

        StatEq.setZero();

        // Radiative transitions
        for (Size k = 0; k < linedata.nrad; k++)
        {
            const Real v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
            const Real v_JI =                 linedata.Ba[k] * Jeff[p][k];

            // Note: we define our transition matrix as the transpose of R in the paper.
            const Size I = linedata.irad[k];
            const Size J = linedata.jrad[k];

            StatEq(J, I) += v_IJ;
            StatEq(J, J) -= v_JI;
            StatEq(I, J) += v_JI;
            StatEq(I, I) -= v_IJ;
        }

        // Approximated Lambda operator
        for (Size k = 0; k < linedata.nrad; k++)
        {
            for (Size m = 0; m < lambda.get_size(p,k); m++)
            {
                const Size   nr =  lambda.get_nr(p, k, m);
                const Real v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                if (nr != p)
                {
                    throw std::runtime_error ("ERROR: non-local Approximated Lambda operator.");
                }

                // Note: we define our transition matrix as the transpose of R in the paper.
                const Size I = linedata.irad[k];
                const Size J = linedata.jrad[k];

                StatEq(J, I) += v_IJ;
                StatEq(I, I) -= v_IJ;
            }
        }

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
                const Size I = colpar.icol[k];
                const Size J = colpar.jcol[k];

                StatEq(J, I) += v_IJ;
                StatEq(J, J) -= v_JI;
                StatEq(I, J) += v_JI;
                StatEq(I, I) -= v_IJ;
            }
        }

        // Replace the last row with normalisation
        StatEq.row(linedata.nlev-1).setOnes();

        y[linedata.nlev-1] = population_tot[p];


        // Solve statistical equilibrium and store level populations
        population.segment(p*linedata.nlev, linedata.nlev) = StatEq.colPivHouseholderQr().solve(y);

    }) // for all cells

}


inline void LineProducingSpecies :: PORTAL_solve_statistical_equilibrium (
    const Double2      &abundance,
    const Vector<Real> &temperature )
{

    rho.resize(parameters->npoints());

    for (Size p = 0; p < parameters->npoints(); p++)
    {
        rho[p].resize(nalign);
    }


    pc::multi_threading::ThreadPrivate<MatrixXd> mtrx_;
    pc::multi_threading::ThreadPrivate<VectorXd>    y_;
    pc::multi_threading::ThreadPrivate<VectorXd>    x_;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
       mtrx_(i).resize(nalign, nalign);

       y_(i)    = VectorXd::Zero(nalign);
       y_(i)[0] = 1.0;

       x_(i).resize(nalign);
    }


    threaded_for (p, parameters->npoints(),
    {
        MatrixXd& mtrx = mtrx_();
        VectorXd&    y =    y_();
        VectorXd&    x =    x_();

        mtrx.setZero();

        // BL:
        //     adding the radiative transitions to the SEE
        //     matrix.
        for (Size t = 0; t < linedata.nrad; t++)
        {
            const Size k = linedata.irad[t];
            const Size l = linedata.jrad[t];

            // BL: we require the angular momentum states.
            //     jk,jl = 2*angular momentum, so that also half-integer
            //     angular momenta can be treated.
            const Size jk = j_lev[k];
            const Size jl = j_lev[l];

            // BL: if j < 1, then no alignment state (k=2) is possible.
            Size km1 = 1;
            Size km2 = 1;

            if (jk > 1) {km1 = 2;}
            if (jl > 1) {km2 = 2;}

            // BL: transitions j1 > j2
            //                (jk) (jl)
            for (Size k1 = 0; k1 < km1; k1++)
            {
                const Size k_a = a_lev[k][k1];

                // BL: spontaneous emission jk --> jl
                mtrx(k_a,k_a) -= linedata.A[t];

                for (Size k2 = 0; k2 < km1; k2++)
                {
                    const Size k_a2 = a_lev[k][k2];

                    // BL: stimulated emission jk --> jl
                    mtrx(k_a,k_a2) -= rm[0][t][k1][k2][0] * J0[p][t];
                    mtrx(k_a,k_a2) -= rm[0][t][k1][k2][1] * J2[p][t];
                }

                for (Size k2 = 0; k2 < km2; k2++)
                {
                    const Size k_l2 = a_lev[l][k2];

                    // BL: stimulated absorption jl --> jk
                    mtrx(k_a,k_l2) += rp[0][t][k1][k2][0] * J0[p][t];
                    mtrx(k_a,k_l2) += rp[0][t][k1][k2][1] * J2[p][t];
                }
            }

            // BL: transitions j2 > j1
            //                (jk) (jl)
            for (Size k1 = 0; k1 < km2; k1++)
            {
                const Size k_l = a_lev[l][k1];

                for (Size k2 = 0; k2 < km2; k2++)
                {
                    const Size k_l2 = a_lev[l][k2];

                    // BL: stimulated absorption jk --> jl
                    mtrx(k_l,k_l2) -= rm[1][t][k1][k2][0] * J0[p][t];
                    mtrx(k_l,k_l2) -= rm[1][t][k1][k2][1] * J2[p][t];
                }

                for (Size k2 = 0; k2 < km1; k2++)
                {
                    const Size k_a2 = a_lev[k][k2];

                    // BL: stimulated emission jl --> jk
                    mtrx(k_l,k_a2) += rp[1][t][k1][k2][0] * J0[p][t];
                    mtrx(k_l,k_a2) += rp[1][t][k1][k2][1] * J2[p][t];
                }

                if ((km1 >= km2) || (k1 == 0))
                {
                    const Size k_a = a_lev[k][k1];

                    // BL: spontaneous emission jl --> jk
                    mtrx(k_l,k_a) += tp[t][k1];
                }
            }
        }


        // BL:
        //     We add the collisional contributions to the rate matrix
        //     matrices while being sensitive to the alignment decay
        //     of states for all K's. Off-diagonal elements only for
        //     the K=0, due to not accoutning for C_K, K>2 collisional
        //     elements. See Eq.(7.101) L&L04.

        for (CollisionPartner &colpar : linedata.colpar)
        {
            Real abn = abundance[p][colpar.num_col_partner];
            Real tmp = temperature[p];

            colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
            colpar.interpolate_collision_coefficients (tmp);

            for (Size t = 0; t < colpar.ncol; t++)
            {
                const Size k  = colpar.icol[t];
                const Size l  = colpar.jcol[t];

                const Size jk = j_lev[k];
                const Size jl = j_lev[l];

                // BL: inelastic jl --> jk
                // BL: superelastic jk --> jl

                Size k_a = a_lev[k][0];
                Size k_l = a_lev[l][0];

                const Real g = std::sqrt((jl+1.0) / (1.0+jk));

                mtrx(k_a,k_l) += g * colpar.Ce_intpld[t] * abn;
                mtrx(k_l,k_a) +=     colpar.Cd_intpld[t] * abn / g;

                // BL: if j < 1, then no alignment state (k=2) is possible.
                Size km1 = 1;
                Size km2 = 1;

                if (jk > 1) {km1 = 2;}
                if (jl > 1) {km2 = 2;}

                // BL: transitions j1 > j2
                //                (jk) (jl)
                for (Size k1 = 0; k1 < km1; k1++)
                {
                    k_a = a_lev[k][k1];
                    // BL: superelastic jk --> jl
                    mtrx(k_a,k_a) -= colpar.Cd_intpld[t] * abn;
                }

                // BL: transitions j2 > j1
                //                (jk) (jl)
                for (Size k1 = 0; k1 < km2; k1++)
                {
                    k_l = a_lev[l][k1];
                    // BL: inelastic jl --> jk
                    mtrx(k_l,k_l) -= colpar.Ce_intpld[t] * abn;
                }
            }
        }

        // BL:
        //     Normalization is added on the first line.
        //     Instead of normalizing using the physical constraint:
        //     \sum_i n_i = 1,
        //     which is usual in the isotropic SEE, we instead normalize the rank-0 irreducible elements
        //     p_i^0 = n_i / sqrt(g_i),
        //     so that
        //     \sum_i sqrt(g_i) p_i = 1
        for (Size s = 0; s < nalign; s++)
        {
            mtrx(0,s) = 0.0;
        }

        for (Size s = 0; s < linedata.nlev; s++)
        {
            const Size k_a = a_lev[s][0];

            mtrx(0,k_a) = std::sqrt(j_lev[s] + 1.0);
        }


        // Solve statistical equilibrium and store level populations
        x = mtrx.colPivHouseholderQr().solve(y);

        for (Size s = 0; s < nalign; s++)
        {
            rho[p][s] = x[s];
        }

    }) // for all points
}


inline void LineProducingSpecies :: PORTAL_solve_statistical_equilibrium_for_point (
    const Double2      &abundance,
    const Vector<Real> &temperature,
    const Size          p           )
{

    rho.resize(parameters->npoints());

    for (Size p = 0; p < parameters->npoints(); p++)
    {
        rho[p].resize(nalign);
    }


    pc::multi_threading::ThreadPrivate<MatrixXd> mtrx_;
    pc::multi_threading::ThreadPrivate<VectorXd>    y_;
    pc::multi_threading::ThreadPrivate<VectorXd>    x_;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
       mtrx_(i).resize(nalign, nalign);

       y_(i)    = VectorXd::Zero(nalign);
       y_(i)[0] = 1.0;

       x_(i).resize(nalign);
    }


    // threaded_for (p, parameters->npoints(),
    // {
        MatrixXd& mtrx = mtrx_();
        VectorXd&    y =    y_();
        VectorXd&    x =    x_();

        mtrx.setZero();

        // BL:
        //     adding the radiative transitions to the SEE
        //     matrix.
        for (Size t = 0; t < linedata.nrad; t++)
        {
            const Size k = linedata.irad[t];
            const Size l = linedata.jrad[t];

            // BL: we require the angular momentum states.
            //     jk,jl = 2*angular momentum, so that also half-integer
            //     angular momenta can be treated.
            const Size jk = j_lev[k];
            const Size jl = j_lev[l];

            // BL: if j < 1, then no alignment state (k=2) is possible.
            Size km1 = 1;
            Size km2 = 1;

            if (jk > 1) {km1 = 2;}
            if (jl > 1) {km2 = 2;}

            // BL: transitions j1 > j2
            //                (jk) (jl)
            for (Size k1 = 0; k1 < km1; k1++)
            {
                const Size k_a = a_lev[k][k1];

                // BL: spontaneous emission jk --> jl
                mtrx(k_a,k_a) -= linedata.A[t];

                for (Size k2 = 0; k2 < km1; k2++)
                {
                    const Size k_a2 = a_lev[k][k2];

                    // BL: stimulated emission jk --> jl
                    mtrx(k_a,k_a2) -= rm[0][t][k1][k2][0] * J0[p][t];
                    mtrx(k_a,k_a2) -= rm[0][t][k1][k2][1] * J2[p][t];
                }

                for (Size k2 = 0; k2 < km2; k2++)
                {
                    const Size k_l2 = a_lev[l][k2];

                    // BL: stimulated absorption jl --> jk
                    mtrx(k_a,k_l2) += rp[0][t][k1][k2][0] * J0[p][t];
                    mtrx(k_a,k_l2) += rp[0][t][k1][k2][1] * J2[p][t];
                }
            }

            // BL: transitions j2 > j1
            //                (jk) (jl)
            for (Size k1 = 0; k1 < km2; k1++)
            {
                const Size k_l = a_lev[l][k1];

                for (Size k2 = 0; k2 < km2; k2++)
                {
                    const Size k_l2 = a_lev[l][k2];

                    // BL: stimulated absorption jk --> jl
                    mtrx(k_l,k_l2) -= rm[1][t][k1][k2][0] * J0[p][t];
                    mtrx(k_l,k_l2) -= rm[1][t][k1][k2][1] * J2[p][t];
                }

                for (Size k2 = 0; k2 < km1; k2++)
                {
                    const Size k_a2 = a_lev[k][k2];

                    // BL: stimulated emission jl --> jk
                    mtrx(k_l,k_a2) += rp[1][t][k1][k2][0] * J0[p][t];
                    mtrx(k_l,k_a2) += rp[1][t][k1][k2][1] * J2[p][t];
                }

                if ((km1 >= km2) || (k1 == 0))
                {
                    const Size k_a = a_lev[k][k1];

                    // BL: spontaneous emission jl --> jk
                    mtrx(k_l,k_a) += tp[t][k1];
                }
            }
        }


        // BL:
        //     We add the collisional contributions to the rate matrix
        //     matrices while being sensitive to the alignment decay
        //     of states for all K's. Off-diagonal elements only for
        //     the K=0, due to not accoutning for C_K, K>2 collisional
        //     elements. See Eq.(7.101) L&L04.

        for (CollisionPartner &colpar : linedata.colpar)
        {
            Real abn = abundance[p][colpar.num_col_partner];
            Real tmp = temperature[p];

            colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
            colpar.interpolate_collision_coefficients (tmp);

            cout << "abn = " << abn << endl;

            for (Size t = 0; t < colpar.ncol; t++)
            {
                const Size k  = colpar.icol[t];
                const Size l  = colpar.jcol[t];

                const Size jk = j_lev[k];
                const Size jl = j_lev[l];

                // BL: inelastic jl --> jk
                // BL: superelastic jk --> jl

                Size k_a = a_lev[k][0];
                Size k_l = a_lev[l][0];

                const Real g = std::sqrt((jl+1.0) / (1.0+jk));

                mtrx(k_a,k_l) += g * colpar.Ce_intpld[t] * abn;
                mtrx(k_l,k_a) +=     colpar.Cd_intpld[t] * abn / g;

                cout << k_a << " " << k_l << " " << g * colpar.Ce_intpld[t] * abn     << endl;
                cout << k_l << " " << k_a << " " <<     colpar.Cd_intpld[t] * abn / g << endl;

                // BL: if j < 1, then no alignment state (k=2) is possible.
                Size km1 = 1;
                Size km2 = 1;

                if (jk > 1) {km1 = 2;}
                if (jl > 1) {km2 = 2;}

                // BL: transitions j1 > j2
                //                (jk) (jl)
                for (Size k1 = 0; k1 < km1; k1++)
                {
                    k_a = a_lev[k][k1];
                    // BL: superelastic jk --> jl
                    mtrx(k_a,k_a) -= colpar.Cd_intpld[t] * abn;

                    cout << k_a << " " << k_a << " " << colpar.Cd_intpld[t] * abn << endl;
                }

                // BL: transitions j2 > j1
                //                (jk) (jl)
                for (Size k1 = 0; k1 < km2; k1++)
                {
                    k_l = a_lev[l][k1];
                    // BL: inelastic jl --> jk
                    mtrx(k_l,k_l) -= colpar.Ce_intpld[t] * abn;

                    cout << k_l << " " << k_l << " " << colpar.Ce_intpld[t] * abn << endl;
                }
            }
        }

        // BL:
        //     Normalization is added on the first line.
        //     Instead of normalizing using the physical constraint:
        //     \sum_i n_i = 1,
        //     which is usual in the isotropic SEE, we instead normalize the rank-0 irreducible elements
        //     p_i^0 = n_i / sqrt(g_i),
        //     so that
        //     \sum_i sqrt(g_i) p_i = 1
        for (Size s = 0; s < nalign; s++)
        {
            mtrx(0,s) = 0.0;
        }

        for (Size s = 0; s < linedata.nlev; s++)
        {
            const Size k_a = a_lev[s][0];

            mtrx(0,k_a) = std::sqrt(j_lev[s] + 1.0);
        }


        // Solve statistical equilibrium and store level populations
        x = mtrx.colPivHouseholderQr().solve(y);

        for (Size s = 0; s < nalign; s++)
        {
            rho[p][s] = x[s];
        }


        M = mtrx;

    // }) // for all points
}
