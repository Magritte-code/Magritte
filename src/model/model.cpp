#include "model.hpp"

#include "paracabs.hpp"
#include "solver/solver.hpp"
#include "tools/heapsort.hpp"

void Model ::read(const Io& io) {
    // Before reading a model, first check whether a file exists at the given
    // location. Otherwise, magritte will try to read nonexistent information and
    // thus segfault.
    if (!io.file_exists(parameters->model_name())) {
        std::cout << "Model file cannot be read." << std::endl;
        std::cout << "No file exists with path: " << parameters->model_name() << std::endl;
        std::cout << "Please check whether the path is spelled correctly." << std::endl;
        throw std::ios_base::failure("Model cannot be found at the specified location");
    }

    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Reading Model...                         " << endl;
    cout << "-------------------------------------------" << endl;
    cout << " model file = " << io.io_file << endl;
    cout << "-------------------------------------------" << endl;

    parameters->read(io);
    geometry.read(io);
    chemistry.read(io);
    thermodynamics.read(io);
    lines.read(io);
    radiation.read(io);

    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Model read, parameters:                  " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  npoints    = " << parameters->npoints() << endl;
    cout << "  nrays      = " << parameters->nrays() << endl;
    cout << "  nboundary  = " << parameters->nboundary() << endl;
    cout << "  nfreqs     = " << parameters->nfreqs() << endl;
    cout << "  nspecs     = " << parameters->nspecs() << endl;
    cout << "  nlspecs    = " << parameters->nlspecs() << endl;
    cout << "  nlines     = " << parameters->nlines() << endl;
    cout << "  nquads     = " << parameters->nquads() << endl;
    cout << "-------------------------------------------" << endl;
    cout << "                                           " << endl;
}

void Model ::write(const Io& io) const {
    // Let only root (rank 0) process write output
    if (pc::message_passing::comm_rank() == 0) {
        parameters->write(io);
        geometry.write(io);
        chemistry.write(io);
        thermodynamics.write(io);
        lines.write(io);
        radiation.write(io);
    }
}

int Model ::compute_inverse_line_widths() {
    cout << "Computing inverse line widths..." << endl;

    lines.set_inverse_width(thermodynamics);

    return (0);
}

///  Compute spectral (=frequency) discretisation
/////////////////////////////////////////////////
int Model ::compute_spectral_discretisation() {
    cout << "Computing spectral discretisation..." << endl;
    radiation.frequencies.resize_data(parameters->nlines() * parameters->nquads());

    threaded_for(p, parameters->npoints(), {
        Real1 freqs(parameters->nfreqs());
        // Size1 nmbrs (parameters->nfreqs());

        Size index0 = 0;
        Size index1 = 0;

        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters->nlspecs(); l++) {
            const Real inverse_mass = lines.lineProducingSpecies[l].linedata.inverse_mass;

            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++) {
                const Real freqs_line = lines.line[index0];
                const Real width      = freqs_line * thermodynamics.profile_width(inverse_mass, p);

                for (Size z = 0; z < parameters->nquads(); z++) {
                    const Real root = lines.lineProducingSpecies[l].quadrature.roots[z];

                    freqs[index1] = freqs_line + width * root;
                    // nmbrs[index1] = index1;

                    index1++;
                }

                index0++;
            }
        }

        // Sort frequencies
        // heapsort (freqs, nmbrs);
        // Warning: when sorting the frequencies, no absolute ordering of the
        // frequency in blocks of the individual lines can be possible
        //  if one encounters overlapping lines. This will ruin some helper
        //  statistics, so it will be disabled from now on.

        // Set all frequencies nu
        for (Size fl = 0; fl < parameters->nfreqs(); fl++) {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }

        // Create lookup table for the frequency corresponding to each line
        // Size1 nmbrs_inverted (parameters->nfreqs());

        for (Size fl = 0; fl < parameters->nfreqs(); fl++) {
            // nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;
            ;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters->nfreqs();
        }

        Size index2 = 0;

        for (Size l = 0; l < parameters->nlspecs(); l++) {
            for (Size k = 0; k < lines.lineProducingSpecies[l].nr_line[p].size(); k++) {
                for (Size z = 0; z < lines.lineProducingSpecies[l].nr_line[p][k].size(); z++) {
                    // NOTE: as the frequencies are no longer sorted, we do not need
                    // to include the point p as index
                    //  lines.lineProducingSpecies[l].nr_line[p][k][z] =
                    //  nmbrs_inverted[index2];
                    lines.lineProducingSpecies[l].nr_line[p][k][z] = index2;

                    radiation.frequencies.appears_in_line_integral[index2] = true;
                    radiation.frequencies.corresponding_l_for_spec[index2] = l;
                    radiation.frequencies.corresponding_k_for_tran[index2] = k;
                    radiation.frequencies.corresponding_z_for_line[index2] = z;

                    index2++;
                }
            }
        }
    })

    // Set corresponding frequencies
    for (Size f = 0; f < parameters->nfreqs(); f++) {
        const Size l = radiation.frequencies.corresponding_l_for_spec[f];
        const Size k = radiation.frequencies.corresponding_k_for_tran[f];

        radiation.frequencies.corresponding_line[f] = lines.line_index(l, k);
    }

    // Set spectral discretisation setting
    spectralDiscretisation = SD_Lines;

    return (0);
}

///  Computer for spectral (=frequency) discretisation
///  Gives same frequency bins to each point
///    @param[in] width : corresponding line width for frequency bins
/////////////////////////////////////////////////////////////////////
int Model ::compute_spectral_discretisation(const Real width) {
    cout << "Computing spectral discretisation..." << endl;
    radiation.frequencies.resize_data(parameters->nlines() * parameters->nquads());

    threaded_for(p, parameters->npoints(), {
        Real1 freqs(parameters->nfreqs());
        // Size1 nmbrs (parameters->nfreqs());

        Size index0 = 0;
        Size index1 = 0;

        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters->nlspecs(); l++) {
            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++) {
                const Real freqs_line = lines.line[index0];

                for (Size z = 0; z < parameters->nquads(); z++) {
                    const Real root = lines.lineProducingSpecies[l].quadrature.roots[z];

                    freqs[index1] = freqs_line + freqs_line * width * root;
                    // nmbrs[index1] = index1;

                    index1++;
                }

                index0++;
            }
        }

        // Sort frequencies
        // heapsort (freqs, nmbrs);

        // Set all frequencies nu
        for (Size fl = 0; fl < parameters->nfreqs(); fl++) {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }

        // Create lookup table for the frequency corresponding to each line
        // Size1 nmbrs_inverted (parameters->nfreqs());

        for (Size fl = 0; fl < parameters->nfreqs(); fl++) {
            // nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;
            ;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters->nfreqs();
        }

        Size index2 = 0;

        for (Size l = 0; l < parameters->nlspecs(); l++) {
            for (Size k = 0; k < lines.lineProducingSpecies[l].nr_line[p].size(); k++) {
                for (Size z = 0; z < lines.lineProducingSpecies[l].nr_line[p][k].size(); z++) {
                    // lines.lineProducingSpecies[l].nr_line[p][k][z] =
                    // nmbrs_inverted[index2];
                    lines.lineProducingSpecies[l].nr_line[p][k][z] = index2;

                    radiation.frequencies.appears_in_line_integral[index2] = true;
                    radiation.frequencies.corresponding_l_for_spec[index2] = l;
                    radiation.frequencies.corresponding_k_for_tran[index2] = k;
                    radiation.frequencies.corresponding_z_for_line[index2] = z;

                    index2++;
                }
            }
        }
    })

    // Set spectral discretisation setting
    spectralDiscretisation = SD_Image;

    return (0);
}

///  Wrapper for compute_spectral_discretisation, filling in
///  parameters.nlines()*parameters.nquads() as number of frequencies, similar
///  to old behavior
int Model ::compute_spectral_discretisation(const Real nu_min, const Real nu_max) {
    return compute_spectral_discretisation(
        nu_min, nu_max, parameters->nlines() * parameters->nquads());
}

///  Computer for spectral (=frequency) discretisation
///  Gives same frequency bins to each point
///    @param[in] min : minimal frequency
///    @param[in] max : maximal frequency
///    @param[in] n_image_freqs : number of frequencies in the discretization
///////////////////////////////////////////////////////
int Model ::compute_spectral_discretisation(const Real nu_min, const Real nu_max,
    const Size n_image_freqs) // TODO: or use nquads() instead?
{
    if (n_image_freqs < 1) {
        throw std::runtime_error("At least a single frequency is needed to compute "
                                 "a spectral discretization.");
    }

    radiation.frequencies.resize_data(n_image_freqs);
    cout << "Computing spectral discretisation..." << endl;

    if (n_image_freqs == 1) { // avoiding division by 0
        threaded_for(p, parameters->npoints(), {
            radiation.frequencies.nu(p, 0) = (Real)(nu_min + nu_max) / 2.0;

            radiation.frequencies.appears_in_line_integral[0] = false;
            radiation.frequencies.corresponding_l_for_spec[0] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[0] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[0] = parameters->nfreqs();
        })
    } else {
        const long double dnu = (nu_max - nu_min) / (n_image_freqs - 1);

        threaded_for(p, parameters->npoints(), {
            for (Size f = 0; f < n_image_freqs; f++) {
                radiation.frequencies.nu(p, f) = (Real)(nu_min + f * dnu);

                radiation.frequencies.appears_in_line_integral[f] = false;
                radiation.frequencies.corresponding_l_for_spec[f] = parameters->nfreqs();
                radiation.frequencies.corresponding_k_for_tran[f] = parameters->nfreqs();
                radiation.frequencies.corresponding_z_for_line[f] = parameters->nfreqs();
            }
        })
    }
    // Set spectral discretisation setting
    spectralDiscretisation = SD_Image;

    // TODO: for all frequencies, set the correct corresponding line. In this way,
    // the OneLine approximation will be compatible with the imagers

    return (0);
}

/// Computer for the level populations, assuming LTE
////////////////////////////////////////////////////
int Model ::compute_LTE_level_populations() {
    cout << "Computing LTE level populations..." << endl;

    // Initialize levels, emissivities and opacities with LTE values
    lines.iteration_using_LTE(chemistry.species.abundance, thermodynamics.temperature.gas);

    return (0);
}

///  Computer for the radiation field
/////////////////////////////////////
int Model ::compute_radiation_field_shortchar_order_0() {
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup<CoMoving>(*this);

    if (parameters->one_line_approximation) {
        solver.solve_shortchar_order_0<OneLine>(*this);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.solve_shortchar_order_0<None>(*this);
        return (0);
    }

    solver.solve_shortchar_order_0<CloseLines>(*this);
    return (0);
}

///  Computer for the radiation field
/////////////////////////////////////
int Model ::compute_radiation_field_feautrier_order_2() {
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup<CoMoving>(*this);

    if (parameters->one_line_approximation) {
        solver.solve_feautrier_order_2<OneLine>(*this);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.solve_feautrier_order_2<None>(*this);
        return (0);
    }

    solver.solve_feautrier_order_2<CloseLines>(*this);
    return (0);
}

/// BUGGED: v computation is incorrect
// ///  Computer for the radiation field
// /////////////////////////////////////
// int Model :: compute_radiation_field_feautrier_order_2_uv ()
// {
//     cout << "Computing radiation field..." << endl;
//
//     Solver solver;
//     solver.setup <CoMoving>                  (*this);
//
//     if (parameters->one_line_approximation)
//     {
//         solver.solve_feautrier_order_2_uv <OneLine> (*this);
//         return (0);
//     }
//
//     if (parameters->sum_opacity_emissivity_over_all_lines)
//     {
//         solver.solve_feautrier_order_2_uv <None> (*this);
//         return (0);
//     }
//
//     solver.solve_feautrier_order_2_uv <CloseLines> (*this);
//     return (0);
// }

///  Computer for the radiation field
/////////////////////////////////////
int Model ::compute_radiation_field_feautrier_order_2_anis() {
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup<CoMoving>(*this);

    if (parameters->one_line_approximation) {
        solver.solve_feautrier_order_2_anis<OneLine>(*this);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.solve_feautrier_order_2_anis<None>(*this);
        return (0);
    }

    solver.solve_feautrier_order_2_anis<CloseLines>(*this);
    return (0);
}

///  Computer for the radiation field
/////////////////////////////////////
int Model ::compute_radiation_field_feautrier_order_2_sparse() {
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup<CoMoving>(*this);

    if (parameters->one_line_approximation) {
        solver.solve_feautrier_order_2_sparse<OneLine>(*this);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.solve_feautrier_order_2_sparse<None>(*this);
        return (0);
    }

    solver.solve_feautrier_order_2_sparse<CloseLines>(*this);
    return (0);
}

///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model ::compute_Jeff() {
    for (LineProducingSpecies& lspec : lines.lineProducingSpecies) {
        threaded_for(p, parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                const Size1 freq_nrs = lspec.nr_line[p][k];

                // Initialize values
                lspec.Jlin[p][k] = 0.0;

                // Integrate over the line
                for (Size z = 0; z < parameters->nquads(); z++) {
                    lspec.Jlin[p][k] += lspec.quadrature.weights[z] * radiation.J(p, freq_nrs[z]);
                }

                double diff = 0.0;

                // Collect the approximated part
                for (Size m = 0; m < lspec.lambda.get_size(p, k); m++) {
                    const Size I =
                        lspec.index(lspec.lambda.get_nr(p, k, m), lspec.linedata.irad[k]);

                    diff += lspec.lambda.get_Ls(p, k, m) * lspec.population[I];
                }

                lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
                lspec.Jdif[p][k] = HH_OVER_FOUR_PI * diff;
            }
        })
    }

    return (0);
}

///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model ::compute_Jeff_sparse() {
    for (LineProducingSpecies& lspec : lines.lineProducingSpecies) {
        threaded_for(p, parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                double diff = 0.0;

                // Collect the approximated part
                for (Size m = 0; m < lspec.lambda.get_size(p, k); m++) {
                    const Size I =
                        lspec.index(lspec.lambda.get_nr(p, k, m), lspec.linedata.irad[k]);

                    diff += lspec.lambda.get_Ls(p, k, m) * lspec.population[I];
                }

                lspec.Jlin[p][k] = lspec.J(p, k);
                lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
            }
        })
    }

    return (0);
}

///  compute level populations from statistical equilibrium
///////////////////////////////////////////////////////////
int Model ::compute_level_populations_from_stateq() {
    lines.iteration_using_statistical_equilibrium(
        chemistry.species.abundance, thermodynamics.temperature.gas, parameters->pop_prec);

    return (0);
}

// Default ng-acceleration used after every
// 'parameters->Ng_acceleration_mem_limit' normal iteration steps
template <>
std::tuple<bool, Size> Model ::ng_acceleration_criterion<Default>(
    bool use_Ng_acceleration, Size prior_normal_iterations) {
    // Ng acceleration will not be used if the user does not allow it
    if (use_Ng_acceleration == false) {
        return std::make_tuple(false, 0);
    }

    const Size skip_N_its = parameters->Ng_acceleration_remove_N_its;
    if (prior_normal_iterations == (parameters->Ng_acceleration_mem_limit + skip_N_its)) {
        return std::make_tuple(true, parameters->Ng_acceleration_mem_limit);
    }
    return std::make_tuple(false, 0);
}

// Adaptive ng acceleration, similar to as described in F. De Ceusters thesis
template <>
std::tuple<bool, Size> Model ::ng_acceleration_criterion<Adaptive>(
    bool use_Ng_acceleration, Size prior_normal_iterations) {
    const Size skip_N_its = parameters->Ng_acceleration_remove_N_its;
    // Ng acceleration will not be used if we have no prior data (or the user does
    // not allow it) For useful acceleration steps, the order needs to be at least
    // 2
    if (use_Ng_acceleration == false
        || prior_normal_iterations < parameters->adaptive_Ng_acceleration_min_order + skip_N_its) {
        return std::make_tuple(false, 0);
    }

    // Computing the current relative change to decide whether or not to
    // use an early Ng-acceleration step
    double sum_fnc_curr = 0.0;
    for (int l = 0; l < parameters->nlspecs(); l++) {
        const double fnc = parameters->adaptive_Ng_acceleration_use_max_criterion
                             ? lines.lineProducingSpecies[l].relative_change_max
                             : lines.lineProducingSpecies[l].relative_change_mean;
        sum_fnc_curr += fnc;
    }

    // the number of previous iterations to use is bounded by the memory limit
    const Size nb_prev_iterations = prior_normal_iterations;

    lines.trial_iteration_using_adaptive_Ng_acceleration(
        parameters->pop_prec, nb_prev_iterations - skip_N_its);

    // Computing the Ng-accelerated relative change to decide whether or
    // not to use an early Ng-acceleration step; compares against previous acceleration step
    double sum_fnc_ng = 0.0;
    for (int l = 0; l < parameters->nlspecs(); l++) {
        const double fnc = parameters->adaptive_Ng_acceleration_use_max_criterion
                             ? lines.lineProducingSpecies[l].relative_change_max
                             : lines.lineProducingSpecies[l].relative_change_mean;
        sum_fnc_ng += fnc;
    }

    // If the acceleration steps start to converge slower than default iterations, it is to to
    // finally apply ng-accelerations
    if ((sum_fnc_ng < sum_fnc_curr)
        || prior_normal_iterations == (parameters->Ng_acceleration_mem_limit + skip_N_its)) {
        return std::make_tuple(true, nb_prev_iterations - skip_N_its);
    }

    return std::make_tuple(false, 0);
}

///  Compute level populations self-consistenly with the radiation field
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] io                  : io object (for writing level populations)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return number of iteration done
///////////////////////////////////////////////////////////////////////////////
int Model ::compute_level_populations(const bool use_Ng_acceleration, const long max_niterations) {
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines) {
        throw std::runtime_error("Spectral discretisation was not set for Lines!");
    }

    // Initialize the number of iterations
    int iteration        = 0;
    int iteration_normal = 0;

    // Initialize errors
    error_mean.clear();
    error_max.clear();

    // Initialize some_not_converged
    bool some_not_converged = true;

    // Iterate as long as some levels are not converged
    while (some_not_converged && (iteration < max_niterations)) {

        std::tuple<bool, Size> tuple_ng_decision;
        if (parameters->use_adaptive_Ng_acceleration) {
            tuple_ng_decision =
                ng_acceleration_criterion<Adaptive>(use_Ng_acceleration, iteration_normal);
        } else {
            tuple_ng_decision =
                ng_acceleration_criterion<Default>(use_Ng_acceleration, iteration_normal);
        }
        const bool use_ng_acceleration_step = std::get<0>(tuple_ng_decision);
        const Size ng_acceleration_order    = std::get<1>(tuple_ng_decision);

        std::cout << "using ng acceleration? " << use_ng_acceleration_step << std::endl;
        if (use_ng_acceleration_step) {
            std::cout << "max order: " << iteration_normal
                      << " used order: " << ng_acceleration_order << std::endl;
            lines.iteration_using_Ng_acceleration(chemistry.species.abundance,
                thermodynamics.temperature.gas, parameters->pop_prec, ng_acceleration_order);

            iteration_normal = 0;
        } else {
            iteration++; // only counting default iterations, as the ng-accelerated
                         // iterations are orders of magnitude faster.

            // Start assuming convergence
            some_not_converged = false;

            // logger.write ("Starting iteration ", iteration);
            cout << "Starting iteration " << iteration << endl;

            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_feautrier_order_2();
            compute_Jeff();

            timer_1.stop();
            timer_1.print_total();

            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium(
                chemistry.species.abundance, thermodynamics.temperature.gas, parameters->pop_prec);

            timer_2.stop();
            timer_2.print_total();

            iteration_normal++;

            for (int l = 0; l < parameters->nlspecs(); l++) {
                error_mean.push_back(lines.lineProducingSpecies[l].relative_change_mean);
                error_max.push_back(lines.lineProducingSpecies[l].relative_change_max);

                if (lines.lineProducingSpecies[l].fraction_not_converged
                    > 1.0 - parameters->convergence_fraction) {
                    some_not_converged = true;
                }

                const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

                // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
                cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
            }
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;

    return iteration;
}

///  Compute level populations self-consistenly with the radiation field
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] io                  : io object (for writing level populations)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return number of iteration done
///////////////////////////////////////////////////////////////////////////////
int Model ::compute_level_populations_sparse(
    const bool use_Ng_acceleration, const long max_niterations) {
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines) {
        throw std::runtime_error("Spectral discretisation was not set for Lines!");
    }

    // Initialize the number of iterations
    int iteration        = 0;
    int iteration_normal = 0;

    // Initialize errors
    error_mean.clear();
    error_max.clear();

    // Initialize some_not_converged
    bool some_not_converged = true;

    // Iterate as long as some levels are not converged
    while (some_not_converged && (iteration < max_niterations)) {

        std::tuple<bool, Size> tuple_ng_decision;
        if (parameters->use_adaptive_Ng_acceleration) {
            tuple_ng_decision =
                ng_acceleration_criterion<Adaptive>(use_Ng_acceleration, iteration_normal);
        } else {
            tuple_ng_decision =
                ng_acceleration_criterion<Default>(use_Ng_acceleration, iteration_normal);
        }
        const bool use_ng_acceleration_step = std::get<0>(tuple_ng_decision);
        const Size ng_acceleration_order    = std::get<1>(tuple_ng_decision);

        std::cout << "using ng acceleration? " << use_ng_acceleration_step << std::endl;
        if (use_ng_acceleration_step) {
            std::cout << "Ng acceleration max order: " << iteration_normal
                      << " used order: " << ng_acceleration_order << std::endl;
            lines.iteration_using_Ng_acceleration(chemistry.species.abundance,
                thermodynamics.temperature.gas, parameters->pop_prec, ng_acceleration_order);

            iteration_normal = 0;
        } else {
            iteration++; // only counting default iterations, as the ng-accelerated
                         // iterations are orders of magnitude faster.

            // Start assuming convergence
            some_not_converged = false;

            // logger.write ("Starting iteration ", iteration);
            cout << "Starting iteration " << iteration << endl;

            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_feautrier_order_2_sparse();
            compute_Jeff_sparse();

            timer_1.stop();
            timer_1.print_total();

            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium_sparse(
                chemistry.species.abundance, thermodynamics.temperature.gas, parameters->pop_prec);

            timer_2.stop();
            timer_2.print_total();

            iteration_normal++;

            for (int l = 0; l < parameters->nlspecs(); l++) {
                error_mean.push_back(lines.lineProducingSpecies[l].relative_change_mean);
                error_max.push_back(lines.lineProducingSpecies[l].relative_change_max);

                if (lines.lineProducingSpecies[l].fraction_not_converged
                    > 1.0 - parameters->convergence_fraction) {
                    some_not_converged = true;
                }

                const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

                // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
                cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
            }
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;

    return iteration;
}

///  Compute level populations self-consistenly with the radiation field
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] io                  : io object (for writing level populations)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return number of iteration done
///////////////////////////////////////////////////////////////////////////////
int Model ::compute_level_populations_shortchar(
    const bool use_Ng_acceleration, const long max_niterations) {
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines) {
        throw std::runtime_error("Spectral discretisation was not set for Lines!");
    }

    // Initialize the number of iterations
    int iteration        = 0;
    int iteration_normal = 0;

    // Initialize errors
    error_mean.clear();
    error_max.clear();

    // Initialize some_not_converged
    bool some_not_converged = true;

    // Iterate as long as some levels are not converged
    while (some_not_converged && (iteration < max_niterations)) {

        std::tuple<bool, Size> tuple_ng_decision;
        if (parameters->use_adaptive_Ng_acceleration) {
            tuple_ng_decision =
                ng_acceleration_criterion<Adaptive>(use_Ng_acceleration, iteration_normal);
        } else {
            tuple_ng_decision =
                ng_acceleration_criterion<Default>(use_Ng_acceleration, iteration_normal);
        }

        const bool use_ng_acceleration_step = std::get<0>(tuple_ng_decision);
        const Size ng_acceleration_order    = std::get<1>(tuple_ng_decision);

        if (use_ng_acceleration_step) {
            std::cout << "Ng acceleration max order: " << iteration_normal
                      << " used order: " << ng_acceleration_order << std::endl;
            lines.iteration_using_Ng_acceleration(chemistry.species.abundance,
                thermodynamics.temperature.gas, parameters->pop_prec, ng_acceleration_order);

            iteration_normal = 0;
        } else {
            iteration++; // only counting default iterations, as the ng-accelerated
                         // iterations are orders of magnitude faster.

            // Start assuming convergence
            some_not_converged = false;

            // logger.write ("Starting iteration ", iteration);
            cout << "Starting iteration " << iteration << endl;

            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_shortchar_order_0();
            compute_Jeff();

            timer_1.stop();
            timer_1.print_total();

            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium(
                chemistry.species.abundance, thermodynamics.temperature.gas, parameters->pop_prec);

            timer_2.stop();
            timer_2.print_total();

            iteration_normal++;

            for (int l = 0; l < parameters->nlspecs(); l++) {
                error_mean.push_back(lines.lineProducingSpecies[l].relative_change_mean);
                error_max.push_back(lines.lineProducingSpecies[l].relative_change_max);

                if (lines.lineProducingSpecies[l].fraction_not_converged
                    > 1.0 - parameters->convergence_fraction) {
                    some_not_converged = true;
                }

                const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

                // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
                cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
            }
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;

    return iteration;
}

///  Computer for the radiation field
/////////////////////////////////////
int Model ::compute_image(const Size ray_nr) {
    cout << "Computing image..." << endl;

    Solver solver;
    solver.setup<Rest>(*this);
    if (parameters->one_line_approximation) {
        throw std::runtime_error("One line approximation is currently not supported for imaging.");
        solver.image_feautrier_order_2<OneLine>(*this, ray_nr);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.image_feautrier_order_2<None>(*this, ray_nr);
        return (0);
    }

    solver.image_feautrier_order_2<CloseLines>(*this, ray_nr);
    return (0);
}

// ///  Wrapper for the new imager
// ///////////////////////////////
// int Model :: compute_image_new (const Vector3D raydir)
// {
//     return compute_image_new(raydir, 256, 256);
// }

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_new(const Size ray_nr) {
    return compute_image_new(geometry.rays.direction[ray_nr], 256, 256);
}

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_new(const Size ray_nr, const Size Nxpix, const Size Nypix) {
    return compute_image_new(geometry.rays.direction[ray_nr], Nxpix, Nypix);
}

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_new(
    const double rx, const double ry, const double rz, const Size Nxpix, const Size Nypix) {
    const Vector3D raydir = Vector3D(rx, ry, rz); // will be normed later on (if not yet normed)
    return compute_image_new(raydir, Nxpix, Nypix);
}

///  Computer for the radiation field, using a new imager TODO: check whether
///  direction is correct (I suspect it is not)
/////////////////////////////////////
int Model ::compute_image_new(const Vector3D raydir, const Size Nxpix, const Size Nypix) {
    if (raydir.squaredNorm() == 0.0) {
        throw std::runtime_error(
            "The given ray direction vector does not point in a direction. Please "
            "use a non-zero (normed) direction vector to generate an image.");
    }
    const Vector3D normed_raydir = raydir * (1 / std::sqrt(raydir.squaredNorm()));
    cout << "Computing image new..." << endl;

    Solver solver;
    // setup has to be handled after image creation, due to the rays themselves
    // depend on the image pixels
    //  solver.setup_new_imager <Rest> (*this);//traced ray length might be
    //  different, thus we might need longer data types
    if (parameters->one_line_approximation) {
        throw std::runtime_error("One line approximation is not supported for imaging.");
        solver.image_feautrier_order_2_new_imager<OneLine>(*this, normed_raydir, Nxpix, Nypix);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.image_feautrier_order_2_new_imager<None>(*this, normed_raydir, Nxpix, Nypix);
        return (0);
    }

    solver.image_feautrier_order_2_new_imager<CloseLines>(*this, normed_raydir, Nxpix, Nypix);
    return (0);
}

///  Computer for image in one point
////////////////////////////////////
int Model ::compute_image_for_point(const Size ray_nr, const Size p) {
    cout << "Computing image for point " << p << "..." << endl;

    Solver solver;
    solver.setup<Rest>(*this);
    if (parameters->one_line_approximation) {
        throw std::runtime_error("One line approximation is currently not supported for imaging.");
        solver.image_feautrier_order_2_for_point<OneLine>(*this, ray_nr, p);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.image_feautrier_order_2_for_point<None>(*this, ray_nr, p);
        return (0);
    }

    solver.image_feautrier_order_2_for_point<CloseLines>(*this, ray_nr, p);
    return (0);
}

///  Computer for optical depth image
//////////////////////////////////////
int Model ::compute_image_optical_depth(const Size ray_nr) {
    Solver solver;
    solver.setup<Rest>(*this);
    if (parameters->one_line_approximation) {
        throw std::runtime_error("One line approximation is currently not supported for imaging.");
        solver.image_optical_depth<OneLine>(*this, ray_nr);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.image_optical_depth<None>(*this, ray_nr);
        return (0);
    }

    solver.image_optical_depth<CloseLines>(*this, ray_nr);
    return (0);
}

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_optical_depth_new(const Size ray_nr) {
    return compute_image_optical_depth_new(geometry.rays.direction[ray_nr], 256, 256);
}

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_optical_depth_new(const Size ray_nr, const Size Nxpix, const Size Nypix) {
    return compute_image_optical_depth_new(geometry.rays.direction[ray_nr], Nxpix, Nypix);
}

///  Wrapper for the new imager
///////////////////////////////
int Model ::compute_image_optical_depth_new(
    const double rx, const double ry, const double rz, const Size Nxpix, const Size Nypix) {
    const Vector3D raydir = Vector3D(rx, ry, rz); // will be normed later on (if not yet normed)
    return compute_image_optical_depth_new(raydir, Nxpix, Nypix);
}

///  Computer for new optical depth image
//////////////////////////////////////
int Model ::compute_image_optical_depth_new(
    const Vector3D raydir, const Size Nxpix, const Size Nypix) {
    if (raydir.squaredNorm() == 0.0) {
        throw std::runtime_error(
            "The given ray direction vector does not point in a direction. Please "
            "use a non-zero (normed) direction vector to generate an image.");
    }
    const Vector3D normed_raydir = raydir * (1 / std::sqrt(raydir.squaredNorm()));
    cout << "Computing image optical depth new..." << endl;

    Solver solver;

    if (parameters->one_line_approximation) {
        throw std::runtime_error("One line approximation is currently not supported for imaging.");
        solver.image_optical_depth_new_imager<OneLine>(*this, normed_raydir, Nxpix, Nypix);
        return (0);
    }

    if (parameters->sum_opacity_emissivity_over_all_lines) {
        solver.image_optical_depth_new_imager<None>(*this, normed_raydir, Nxpix, Nypix);
        return (0);
    }

    solver.image_optical_depth_new_imager<CloseLines>(*this, normed_raydir, Nxpix, Nypix);
    return (0);
}

int Model ::set_eta_and_chi(const Size rr) {
    Solver solver;
    solver.set_eta_and_chi(*this, rr);

    return (0);
}

int Model ::set_boundary_condition() {
    Solver solver;
    solver.set_boundary_condition(*this);

    return (0);
}

int Model ::set_column() {
    Solver solver;
    solver.set_column(*this);

    return (0);
}

// /  Setter for the maximum allowed shift value determined by the smallest line
// /////////////////////////////////////////////////////////////////////////////
int Model ::set_dshift_max() {
    // Allocate memory
    dshift_max.resize(parameters->npoints());

    // For all points
    threaded_for(o, parameters->npoints(), {
        dshift_max[o] = std::numeric_limits<Real>::max();

        for (const LineProducingSpecies& lspec : lines.lineProducingSpecies) {
            const Real inverse_mass = lspec.linedata.inverse_mass;
            const Real new_dshift_max =
                parameters->max_width_fraction * thermodynamics.profile_width(inverse_mass, o);

            if (dshift_max[o] > new_dshift_max) {
                dshift_max[o] = new_dshift_max;
            }
        }
    })

    return (0);
}
