#include "paracabs.hpp"
#include "tools/constants.hpp"
#include "tools/types.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>

///  Indexer for level populations
///    @param[in] p : index of the cell
///    @param[in] i : index of the level
///    @return corresponding index for p and i
//////////////////////////////////////////////
inline Size LineProducingSpecies ::index(const Size p, const Size i) const {
    return i + p * linedata.nlev;
}

///  Getter for the line emissivity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line emissivity for cell p and transition k
//////////////////////////////////////////////////////////
inline Real LineProducingSpecies ::get_emissivity(const Size p, const Size k) const {
    const Size i = index(p, linedata.irad[k]);

    return HH_OVER_FOUR_PI * linedata.A[k] * population(i);
}

///  Getter for the line opacity
///    @param[in] p : index of the cell
///    @param[in] k : index of the transition
///    @return line opacity for cell p and transition k
///////////////////////////////////////////////////////
inline Real LineProducingSpecies ::get_opacity(const Size p, const Size k) const {
    const Size i = index(p, linedata.irad[k]);
    const Size j = index(p, linedata.jrad[k]);

    return HH_OVER_FOUR_PI * (population(j) * linedata.Ba[k] - population(i) * linedata.Bs[k]);
}

///  set_LTE_level_populations
///    @param[in] abundance: abundance of line species
///    @param[in] temperature: local gas temperature
///    @param[in] p: number of cell
///    @param[in] l: number of line producing species
///////////////////////////////////////////////////////////
inline void LineProducingSpecies ::update_using_LTE(
    const Double2& abundance, const Vector<Real>& temperature) {
    threaded_for(p, parameters->npoints(), {
        population_tot[p] = abundance[p][linedata.num];

        Real partition_function = 0.0;

        for (Size i = 0; i < linedata.nlev; i++) {
            const Size ind = index(p, i);

            population(ind) = linedata.weight[i] * exp(-linedata.energy[i] / (KB * temperature[p]));

            partition_function += population(ind);
        }

        for (Size i = 0; i < linedata.nlev; i++) {
            const Size ind = index(p, i);

            population(ind) *= population_tot[p] / partition_function;
        }
    })

    populations.push_back(population);
}

///  Set a single line, at a given point to LTE
///    @param[in] temperature: local gas temperature
///    @param[in] p: number of cell
///    @param[in] k: number of line
/// Note: do call this from back to front in the loop over k, as otherwise new population inversions
/// might occur
///////////////////////////////////////////////////////////
inline void LineProducingSpecies ::set_LTE(
    const Vector<Real>& temperature, const Size p, const Size k) {
    // population_tot[p] = abundance[p][linedata.num];
    const Size i = index(p, linedata.irad[k]);
    const Size j = index(p, linedata.jrad[k]);

    Real population_tot_line = population(i) + population(j);
    Real partition_function  = 0.0;

    population(i) = linedata.weight[linedata.irad[k]]
                  * exp(-linedata.energy[linedata.irad[k]] / (KB * temperature[p]));
    population(j) = linedata.weight[linedata.jrad[k]]
                  * exp(-linedata.energy[linedata.jrad[k]] / (KB * temperature[p]));
    partition_function = population(i) + population(j);
    population(i) *= population_tot_line / partition_function;
    population(j) *= population_tot_line / partition_function;
}

///  Set the given levels, at a given point to LTE
///    @param[in] temperature: local gas temperature
///    @param[in] p: number of cell
///    @param[in] level_mask: vector of bools, indicating which levels should be set to LTE
inline void LineProducingSpecies ::set_LTE_specific_levels(
    const Vector<Real>& temperature, const Size p, const vector<char>& level_mask) {
    Real population_tot_levels = 0.0;
    Real partition_function    = 0.0;
    for (Size i = 0; i < linedata.nlev; i++) {
        if (level_mask[i]) {
            const Size ind = index(p, i);
            population_tot_levels += population(ind);
            population(ind) = linedata.weight[i] * exp(-linedata.energy[i] / (KB * temperature[p]));
            partition_function += population(ind);
        }
    }
    for (Size i = 0; i < linedata.nlev; i++) {
        if (level_mask[i]) {
            const Size ind = index(p, i);
            population(ind) *= population_tot_levels / partition_function;
        }
    }
}

/// This function check whether the level populations have converged, using the
/// specified precision
inline void LineProducingSpecies ::check_for_convergence(const Real pop_prec) {
    const Real weight = 1.0 / (parameters->npoints() * linedata.nlev);

    Real fnc   = 0.0;
    Real rcm   = 0.0;
    Real rcmax = 0.0;

#pragma omp parallel for reduction(+ : fnc, rcm) reduction(max : rcmax)
    for (Size p = 0; p < parameters->npoints(); p++) {
        for (Size i = 0; i < linedata.nlev; i++) {
            const Size ind = index(p, i);

            if (population(ind) > parameters->min_rel_pop_for_convergence * population_tot[p]) {
                // Real relative_change = 2.0;
                //
                // relative_change *= std::abs(population (ind) - population_prev1
                // (ind)); relative_change /= (population (ind) + population_prev1
                // (ind));
                // minor computed negative level populations might result in negative
                // values for this
                const Real relative_change = 2.0
                                           * std::abs((population(ind) - population_prev1(ind))
                                                      / (population(ind) + population_prev1(ind)));

                if (relative_change > pop_prec) {
                    fnc += weight;
                }

                rcm += (weight * relative_change);
                rcmax = std::max(relative_change, rcmax);
            }
        }
    }

    fraction_not_converged = fnc;
    relative_change_mean   = rcm;
    relative_change_max    = rcmax;
}

/// This function checks whether the trial level populations have converged,
/// using the specified precision
inline void LineProducingSpecies ::check_for_convergence_trial(const Real pop_prec) {
    const Real weight = 1.0 / (parameters->npoints() * linedata.nlev);

    Real fnc   = 0.0;
    Real rcm   = 0.0;
    Real rcmax = 0.0;

#pragma omp parallel for reduction(+ : fnc, rcm) reduction(max : rcmax)
    for (Size p = 0; p < parameters->npoints(); p++) {
        for (Size i = 0; i < linedata.nlev; i++) {
            const Size ind = index(p, i);

            // if (population(ind) > parameters->min_rel_pop_for_convergence *
            // population_tot[p]) { Real relative_change = 2.0;
            //
            // relative_change *= std::abs(trial_population (ind) - trial_population_prev
            // (ind)); relative_change /= (trial_population (ind) + trial_population_prev
            // (ind)); minor computed negative level populations might result in negative values
            // for this
            const Real relative_change =
                2.0
                * std::abs((trial_population(ind) - trial_population_prev(ind))
                           / (trial_population(ind) + trial_population_prev(ind)));

            if (relative_change > pop_prec) {
                fnc += weight;
            }

            rcm += (weight * relative_change);
            rcmax = std::max(relative_change, rcmax);
            // }
        }
    }

    fraction_not_converged = fnc;
    relative_change_mean   = rcm;
    relative_change_max    = rcmax;
}

///  update_using_Ng_acceleration: perform a Ng accelerated iteration step
///    for level populations. All variable names are based on lecture notes
///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
///////////////////////////////////////////////////////////////////////////
void LineProducingSpecies ::update_using_Ng_acceleration() {
    VectorXld Wt(parameters->npoints() * linedata.nlev);

    VectorXld Q1 = population - 2.0 * population_prev1 + population_prev2;
    VectorXld Q2 = population - population_prev1 - population_prev2 + population_prev3;
    VectorXld Q3 = population - population_prev1;

    // OMP_PARALLEL_FOR (ind, ncells*linedata.nlev)
    //{
    //   if (population (ind) > 0.0)
    //   {
    //     Wt (ind) = Jlin[p][k];
    //   }

    //  else
    //  {
    //    Wt (ind) = 1.0;
    //  }
    //}

    // const double A1 = Q1.dot (Wt.asDiagonal()*Q1);
    // const double A2 = Q1.dot (Wt.asDiagonal()*Q2);
    // const double B2 = Q2.dot (Wt.asDiagonal()*Q2);
    // const double C1 = Q1.dot (Wt.asDiagonal()*Q3);
    // const double C2 = Q2.dot (Wt.asDiagonal()*Q3);

    const long double A1 = Q1.dot(Q1);
    const long double A2 = Q1.dot(Q2);
    const long double B2 = Q2.dot(Q2);
    const long double C1 = Q1.dot(Q3);
    const long double C2 = Q2.dot(Q3);

    const long double B1 = A2;

    const long double denominator = A1 * B2 - A2 * B1;

    if (denominator != 0.0) {
        const VectorXld pop_tmp = population;

        const long double a = (C1 * B2 - C2 * B1) / denominator;
        const long double b = (C2 * A1 - C1 * A2) / denominator;

        population = (1.0 - a - b) * population + a * population_prev1 + b * population_prev2;

        population_prev3 = population_prev2;
        population_prev2 = population_prev1;
        population_prev1 = pop_tmp;
    }
}

///  update_using_acceleration: perform a Ng accelerated iteration step
///    for level populations. All variable names are based on lecture notes
///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
///////////////////////////////////////////////////////////////////////////
void LineProducingSpecies ::update_using_acceleration(const Size order) {
    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    residuals.push_back(population - populations.back());
    populations.push_back(population);

    // inequality due to residuals needing to be computed from populations, but
    // populations may be computed without computing residuals
    if (populations.size() < (residuals.size() + 1)) {
        std::cout << "Inconsistent data sizes for ng acceleration." << std::endl;
        std::cout << "Size populations: " << populations.size()
                  << " Size residuals: " << residuals.size() << std::endl;
        throw std::runtime_error("Error during ng acceleration; inconstent sizes.");
    }
    if (residuals.size() < order) {
        std::cout << "Ng acceleration order: " << order << std::endl;
        std::cout << "Size residuals: " << residuals.size() << std::endl;
        throw std::runtime_error("Error during ng acceleration; order of "
                                 "acceleration greather than data stored.");
    }

    const Size popsize = populations.size();
    const Size ressize = residuals.size();

    MatrixXld RTR(order, order);

    for (Size i = 0; i < order; i++) {
        for (Size j = 0; j < order; j++) {
            RTR(i, j) = residuals[ressize - order + i].dot(residuals[ressize - order + j]);
        }
    }

    VectorXld ones = VectorXld::Constant(order, 1.0);
    VectorXld coef = RTR.colPivHouseholderQr().solve(ones);
    coef /= coef.sum();

    population = VectorXld::Zero(population.size());

    for (Size i = 0; i < order; i++) {
        population += populations[popsize - 1 - i] * coef[order - 1 - i];
    }

    // enforcing memory limit by removing almost all previous information after
    // (general) ng-acceleration Last population must be kept to compute residuals
    residuals.clear();
    populations.erase(populations.begin(), populations.end() - 1);
}

///  update_using_acceleration: perform a Ng accelerated iteration step
///    for level populations. All variable names are based on lecture notes
///    by C.P. Dullemond which are based on Olson, Auer and Buchler (1985).
///  Does not update the level populations, storing the results instead into a
///  seperate array
///////////////////////////////////////////////////////////////////////////
void LineProducingSpecies ::update_using_acceleration_trial(const Size order) {
    residuals.push_back(population - populations.back());
    populations.push_back(population);

    trial_population_prev = trial_population; // stores the previous ng acceleration estimate,
                                              // used for the adaptive ng acceleration

    // inequality due to residuals needing to be computed from populations, but
    // populations may be computed without computing residuals
    if (populations.size() < (residuals.size() + 1)) {
        std::cout << "Inconsistent data sizes for ng acceleration." << std::endl;
        std::cout << "Size populations: " << populations.size()
                  << " Size residuals: " << residuals.size() << std::endl;
        throw std::runtime_error("Error during ng acceleration; inconstent sizes.");
    }
    if (populations.size() < order) {
        std::cout << "Ng acceleration order: " << order << std::endl;
        std::cout << "Size residuals: " << residuals.size() << std::endl;
        throw std::runtime_error("Error during ng acceleration; order of "
                                 "acceleration greather than data stored.");
    }

    const Size popsize = populations.size();
    const Size ressize = residuals.size();

    // MatrixXr RTR(order, order);
    MatrixXld RTR(order, order);

    for (Size i = 0; i < order; i++) {
        for (Size j = 0; j < order; j++) {
            RTR(i, j) = residuals[ressize - order + i].dot(residuals[ressize - order + j]);
        }
    }

    VectorXld ones = VectorXld::Constant(order, 1.0);
    VectorXld coef = RTR.colPivHouseholderQr().solve(ones);
    coef /= coef.sum();

    trial_population = VectorXld::Zero(population.size());

    for (Size i = 0; i < order; i++) {
        trial_population += populations[popsize - 1 - i] * coef[order - 1 - i];
    }

    // as this is a trial, the last elements will now be deleted for consistency
    residuals.pop_back();
    populations.pop_back();
}

///  update_using_statistical_equilibrium: computes level populations by solving
///  the statistical equilibrium equation taking into account the radiation
///  field. This function is not parallellized, as it constructs a single large sparse matrix.
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
/////////////////////////////////////////////////////////////////////////////////
inline void LineProducingSpecies::update_using_statistical_equilibrium(
    const Double2& abundance, const Vector<Real>& temperature) {
    RT.resize(parameters->npoints() * linedata.nlev, parameters->npoints() * linedata.nlev);
    LambdaStar.resize(parameters->npoints() * linedata.nlev, parameters->npoints() * linedata.nlev);
    LambdaTest.resize(parameters->npoints() * linedata.nlev, parameters->npoints() * linedata.nlev);

    const Size non_zeros =
        parameters->npoints() * (linedata.nlev + 6 * linedata.nrad + 4 * linedata.ncol_tot);
    // Store previous iterations
    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    residuals.push_back(population - populations.back());
    populations.push_back(population);

    //    SparseMatrix<double> RT (ncells*linedata.nlev, ncells*linedata.nlev);

    VectorXld y = VectorXld::Zero(parameters->npoints() * linedata.nlev);

    vector<Triplet<long double, Size>> triplets;
    //    vector<Triplet<Real, Size>> triplets_LT;
    //    vector<Triplet<Real, Size>> triplets_LS;

    triplets.reserve(non_zeros);
    //    triplets_LT.reserve (non_zeros);
    //    triplets_LS.reserve (non_zeros);

    for (Size p = 0; p < parameters->npoints();
         p++) // !!! no OMP because push_back is not thread safe !!!
    {
        // Radiative transitions

        for (Size k = 0; k < linedata.nrad; k++) {
            const long double v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
            const long double v_JI = linedata.Ba[k] * Jeff[p][k];

            // const Real t_IJ = linedata.Bs[k] * Jdif[p][k];
            // const Real t_JI = linedata.Ba[k] * Jdif[p][k];

            // Note: we define our transition matrix as the transpose of R in the
            // paper.
            const Size I = index(p, linedata.irad[k]);
            const Size J = index(p, linedata.jrad[k]);

            if (linedata.jrad[k] != linedata.nlev - 1) {
                triplets.push_back(Triplet<long double, Size>(J, I, +v_IJ));
                triplets.push_back(Triplet<long double, Size>(J, J, -v_JI));

                // triplets_LS.push_back (Triplet<Real, Size> (J, I, +t_IJ));
                // triplets_LS.push_back (Triplet<Real, Size> (J, J, -t_JI));
            }

            if (linedata.irad[k] != linedata.nlev - 1) {
                triplets.push_back(Triplet<long double, Size>(I, J, +v_JI));
                triplets.push_back(Triplet<long double, Size>(I, I, -v_IJ));

                // triplets_LS.push_back (Triplet<Real, Size> (I, J, +t_JI));
                // triplets_LS.push_back (Triplet<Real, Size> (I, I, -t_IJ));
            }
        }

        // Approximated Lambda operator

        for (Size k = 0; k < linedata.nrad; k++) {
            for (Size m = 0; m < lambda.get_size(p, k); m++) {
                const Size nr          = lambda.get_nr(p, k, m);
                const long double v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                // Note: we define our transition matrix as the transpose of R in the
                // paper.
                const Size I = index(nr, linedata.irad[k]);
                const Size J = index(p, linedata.jrad[k]);

                if (linedata.jrad[k] != linedata.nlev - 1) {
                    triplets.push_back(Triplet<long double, Size>(J, I, +v_IJ));
                    // triplets_LT.push_back (Triplet<Real, Size> (J, I, +v_IJ));
                }

                if (linedata.irad[k] != linedata.nlev - 1) {
                    triplets.push_back(Triplet<long double, Size>(I, I, -v_IJ));
                    // triplets_LT.push_back (Triplet<Real, Size> (I, I, -v_IJ));
                }
            }
        }

        // Collisional transitions

        for (CollisionPartner& colpar : linedata.colpar) {
            Real abn = abundance[p][colpar.num_col_partner];
            Real tmp = temperature[p];

            colpar.adjust_abundance_for_ortho_or_para(tmp, abn);
            colpar.interpolate_collision_coefficients(tmp);

            for (Size k = 0; k < colpar.ncol; k++) {
                const long double v_IJ = colpar.Cd_intpld()[k] * abn;
                const long double v_JI = colpar.Ce_intpld()[k] * abn;

                // Note: we define our transition matrix as the transpose of R in the
                // paper.
                const Size I = index(p, colpar.icol[k]);
                const Size J = index(p, colpar.jcol[k]);

                if (colpar.jcol[k] != linedata.nlev - 1) {
                    triplets.push_back(Triplet<long double, Size>(J, I, +v_IJ));
                    triplets.push_back(Triplet<long double, Size>(J, J, -v_JI));
                }

                if (colpar.icol[k] != linedata.nlev - 1) {
                    triplets.push_back(Triplet<long double, Size>(I, J, +v_JI));
                    triplets.push_back(Triplet<long double, Size>(I, I, -v_IJ));
                }
            }
        }

        for (Size i = 0; i < linedata.nlev; i++) {
            const Size I = index(p, linedata.nlev - 1);
            const Size J = index(p, i);

            triplets.push_back(Triplet<long double, Size>(I, J, 1.0));
        }

        y[index(p, linedata.nlev - 1)] = population_tot[p];

    } // for all cells

    RT.setFromTriplets(triplets.begin(), triplets.end());
    // LambdaStar.setFromTriplets (triplets_LS.begin(), triplets_LS.end());
    // LambdaTest.setFromTriplets (triplets_LT.begin(), triplets_LT.end());

    // cout << "Compressing RT" << endl;

    // RT.makeCompressed ();

    // Eigen::BiCGSTAB <SparseMatrix<double>> solver;

    // cout << "Try compute" << endl;

    // solver.compute (RT);

    // if (solver.info() != Eigen::Success)
    //{
    //   cout << "Decomposition failed" << endl;
    //   //assert(false);
    // }

    // for (int tel=0; tel<5; tel++)
    //{
    //   //Eigen::Gues x0 = population;

    //  population = solver.solveWithGuess (y, population);
    //  std::cout << "#iterations:     " << solver.iterations() << std::endl;
    //  std::cout << "estimated error: " << solver.error()      << std::endl;
    //}

    // assert (false);

    SparseLU<SparseMatrix<long double>, COLAMDOrdering<int>> solver;
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    cout << "Analyzing system of rate equations..." << endl;

    solver.analyzePattern(RT);

    cout << "Factorizing system of rate equations..." << endl;

    solver.factorize(RT);

    if (solver.info() != Eigen::Success) {
        cout << "Factorization failed with error message:" << endl;
        cout << solver.lastErrorMessage() << endl;

        throw std::runtime_error("Eigen solver ERROR.");
    }

    cout << "Solving rate equations for the level populations..." << endl;

    population = solver.solve(y);

    if (solver.info() != Eigen::Success) {
        cout << "Solving failed with error:" << endl;
        cout << solver.lastErrorMessage() << endl;
        assert(false);
    }

    cout << "Succesfully solved for the level populations!" << endl;

    // OMP_PARALLEL_FOR (p, ncells)
    //{
    //
    //   for (long i = 0; i < linedata.nlev; i++)
    //   {
    //     const long I = index (p, i);

    //    population[I] = population_prev1[I];

    //    //if (population[I] < 1.0E-50)
    //    //{
    //    //  population[I] = 1.0E-50;
    //    //}
    //  }
    //}
}

///  update_using_statistical_equilibrium: computes level populations by solving
///  the statistical equilibrium equation taking into account the local radiation
///  field. Is parallelized over the different positions.
///    @param[in] abundance: chemical abundances of species in the model
///    @param[in] temperature: gas temperature in the model
/////////////////////////////////////////////////////////////////////////////////
inline void LineProducingSpecies::update_using_statistical_equilibrium_sparse(
    const Double2& abundance, const Vector<Real>& temperature) {

    if (parameters->n_off_diag != 0) {
        throw std::runtime_error("parameters->n_off_diag != 0 so cannot use sparse update.");
    }

    // Store previous iterations
    population_prev3 = population_prev2;
    population_prev2 = population_prev1;
    population_prev1 = population;

    residuals.push_back(population - populations.back());
    populations.push_back(population);

    pc::multi_threading::ThreadPrivate<MatrixXld> StatEq;
    pc::multi_threading::ThreadPrivate<VectorXld> y;
    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++) {
        StatEq(i).resize(linedata.nlev, linedata.nlev);
        y(i) = VectorXld::Zero(linedata.nlev);
    }

    threaded_for(p,
        parameters->npoints(), //)
                               // for (Size p = 0; p < parameters->npoints();
                               //      p++) // !!! no OMP because push_back is not thread safe !!!
        {
            StatEq().setZero();

            // Radiative transitions
            for (Size k = 0; k < linedata.nrad; k++) {
                const long double v_IJ = linedata.A[k] + linedata.Bs[k] * Jeff[p][k];
                const long double v_JI = linedata.Ba[k] * Jeff[p][k];

                // Note: we define our transition matrix as the transpose of R in the
                // paper.
                const Size I = linedata.irad[k];
                const Size J = linedata.jrad[k];

                StatEq()(J, I) += v_IJ;
                StatEq()(J, J) -= v_JI;
                StatEq()(I, J) += v_JI;
                StatEq()(I, I) -= v_IJ;
            }

            // Approximated Lambda operator
            for (Size k = 0; k < linedata.nrad; k++) {
                for (Size m = 0; m < lambda.get_size(p, k); m++) {
                    const Size nr          = lambda.get_nr(p, k, m);
                    const long double v_IJ = -lambda.get_Ls(p, k, m) * get_opacity(p, k);

                    if (nr != p) {
                        throw std::runtime_error("ERROR: non-local Approximated Lambda operator.");
                    }

                    // Note: we define our transition matrix as the transpose of R in the
                    // paper.
                    const Size I = linedata.irad[k];
                    const Size J = linedata.jrad[k];

                    StatEq()(J, I) += v_IJ;
                    StatEq()(I, I) -= v_IJ;
                }
            }

            // Collisional transitions
            for (CollisionPartner& colpar : linedata.colpar) {
                Real abn = abundance[p][colpar.num_col_partner];
                Real tmp = temperature[p];

                colpar.adjust_abundance_for_ortho_or_para(tmp, abn);
                colpar.interpolate_collision_coefficients(tmp);

                for (Size k = 0; k < colpar.ncol; k++) {
                    const long double v_IJ = colpar.Cd_intpld()[k] * abn;
                    const long double v_JI = colpar.Ce_intpld()[k] * abn;

                    // Note: we define our transition matrix as the transpose of R in the
                    // paper.
                    const Size I = colpar.icol[k];
                    const Size J = colpar.jcol[k];

                    StatEq()(J, I) += v_IJ;
                    StatEq()(J, J) -= v_JI;
                    StatEq()(I, J) += v_JI;
                    StatEq()(I, I) -= v_IJ;
                }
            }

            // Replace the last row with normalisation
            StatEq().row(linedata.nlev - 1).setOnes();

            y()[linedata.nlev - 1] = population_tot[p];

            // Solve statistical equilibrium and store level populations
            population.segment(p * linedata.nlev, linedata.nlev) =
                StatEq().colPivHouseholderQr().solve(y());
        });
    // for all cells

    // OMP_PARALLEL_FOR (p, ncells)
    //{
    //
    //   for (long i = 0; i < linedata.nlev; i++)
    //   {
    //     const long I = index (p, i);

    //    population[I] = population_prev1[I];

    //    //if (population[I] < 1.0E-50)
    //    //{
    //    //  population[I] = 1.0E-50;
    //    //}
    //  }
    //}
}

///  correct_negative_populations: sets negative values in the level populations to zero and
///  renormalizes the other populations. It also sets masering levels to LTE. This function
///  should therefore be called after each level population, as Magritte does not handle
///  negative line opacities.
inline void LineProducingSpecies::correct_negative_populations(
    const Double2& abundance, const Vector<Real>& temperature) {
    threaded_for(p, parameters->npoints(), {
        if (population_tot[p] > 0.0) {
            Real total_positive_population = 0.0;
            for (Size i = 0; i < linedata.nlev; i++) {
                if (population(p * linedata.nlev + i) < 0.0) {
                    population(p * linedata.nlev + i) = 0.0;
                } else {
                    total_positive_population += population(p * linedata.nlev + i);
                }
            }
            // and renormalize the level populations
            for (Size i = 0; i < linedata.nlev;
                 i++) { // TODO: use more fancy operation to divide entire column at once
                        // Also just floor the entire vector by 0
                population(p * linedata.nlev + i) = population(p * linedata.nlev + i)
                                                  * population_tot[p] / total_positive_population;
            }
        }
    });

    bool population_inversions = false;
    // Also check if some population inversions (negative opacities) are present; set these
    // points to LTE instead; This is due to us using the Feautrier solver, which cannot handle
    // negative optical depths, even if no significant masers exist in the model
    threaded_for(p, parameters->npoints(), {
        // keep track of which levels have been modified, in order to do a final pass later
        // use char, because vector<bool> is not a std::vector
        std::vector<char> modified_level(linedata.nlev, false);

        // keep track of how many levels have been modified, as stopping criterion
        Size next_n_modified_levels = 0;
        Size curr_n_modified_levels = 0;

        // We try to keep fixing the levels until they are finally fully consistent...
        // This do while loop is guaranteed to finish, as the number of modified levels can only
        // increase and is bounded from above by linedata.nlev
        do {
            curr_n_modified_levels = next_n_modified_levels;

            bool population_inversion_at_point = false;

            // start by checking transitions one by one to get an idea of which ones need to be set
            // to LTE
            for (Size k = linedata.nrad - 1; k != static_cast<Size>(-1);
                 k--) { // reversed loop, stops at 0
                Size i = index(p, linedata.irad[k]);
                Size j = index(p, linedata.jrad[k]);
                if (population(j) < linedata.Bs[k] / linedata.Ba[k] * population(i)
                                        * parameters->population_inversion_fraction) {
                    set_LTE(temperature, p, k);
                    modified_level[linedata.irad[k]] = true;
                    modified_level[linedata.jrad[k]] = true;
                    population_inversion_at_point    = true;
                    population_inversions = true; // err, technically this isn't thread safe, but we
                                                  // are only trying to set a flag to true
                }
            }
            // And also check from low to high
            for (Size k = 0; k < linedata.nrad; k++) {
                Size i = index(p, linedata.irad[k]);
                Size j = index(p, linedata.jrad[k]);
                if (population(j) < linedata.Bs[k] / linedata.Ba[k] * population(i)
                                        * parameters->population_inversion_fraction) {
                    set_LTE(temperature, p, k);
                    modified_level[linedata.irad[k]] = true;
                    modified_level[linedata.jrad[k]] = true;
                    population_inversion_at_point    = true;
                    population_inversions = true; // err, technically this isn't thread safe, but we
                                                  // are only trying to set a flag to true
                }
            }

            // and finally set all modified levels simultaneously to LTE
            if (population_inversion_at_point) {
                set_LTE_specific_levels(temperature, p, modified_level);
            }

            // and count how many levels have been modified
            next_n_modified_levels =
                std::accumulate(modified_level.begin(), modified_level.end(), static_cast<Size>(0));
        } while (next_n_modified_levels != curr_n_modified_levels);
        // If no additional levels have been modified, we are done
    });

    if (population_inversions) {
        std::cout << "Minor warning: population inversions detected; Magritte does not handle "
                     "masers, so setting affected populations to LTE."
                  << std::endl;
    }
}
