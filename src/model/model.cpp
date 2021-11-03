#include "paracabs.hpp"
#include "model.hpp"
#include "tools/heapsort.hpp"
#include "tools/timer.hpp"
#include "solver/solver.hpp"


void Model :: read (const Io& io)
{
    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Reading Model...                         " << endl;
    cout << "-------------------------------------------" << endl;
    cout << " model file = " << io.io_file                << endl;
    cout << "-------------------------------------------" << endl;

    parameters    .read (io);
    geometry      .read (io);
    chemistry     .read (io);
    thermodynamics.read (io);
    lines         .read (io);
    radiation     .read (io);

    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Model read, parameters:                  " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  npoints    = " << parameters.npoints    () << endl;
    cout << "  nrays      = " << parameters.nrays      () << endl;
    cout << "  nrays_red  = " << parameters.nrays_red  () << endl;
    cout << "  nboundary  = " << parameters.nboundary  () << endl;
    cout << "  nfreqs     = " << parameters.nfreqs     () << endl;
    // cout << "  nfreqs_red = " << parameters.nfreqs_red () << endl;
    cout << "  nspecs     = " << parameters.nspecs     () << endl;
    cout << "  nlspecs    = " << parameters.nlspecs    () << endl;
    cout << "  nlines     = " << parameters.nlines     () << endl;
    cout << "  nquads     = " << parameters.nquads     () << endl;
    cout << "-------------------------------------------" << endl;
    cout << "                                           " << endl;
}


void Model :: write (const Io& io) const
{
    parameters    .write (io);
    geometry      .write (io);
    chemistry     .write (io);
    thermodynamics.write (io);
    lines         .write (io);
    radiation     .write (io);
}


int Model :: compute_inverse_line_widths ()
{
    cout << "Computing inverse line widths..." << endl;

    lines.set_inverse_width (thermodynamics);

    return (0);
}



///  Compute spectral (=frequency) discretisation
/////////////////////////////////////////////////
int Model :: compute_spectral_discretisation ()
{
    cout << "Computing spectral discretisation..." << endl;

    threaded_for (p, parameters.npoints(),
    {
        Real1 freqs (parameters.nfreqs());
        Size1 nmbrs (parameters.nfreqs());

        Size index0 = 0;
        Size index1 = 0;

        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            const Real inverse_mass = lines.lineProducingSpecies[l].linedata.inverse_mass;

            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real freqs_line = lines.line[index0];
                const Real width      = freqs_line * thermodynamics.profile_width (inverse_mass, p);

                for (Size z = 0; z < parameters.nquads(); z++)
                {
                    const Real root = lines.lineProducingSpecies[l].quadrature.roots[z];

                    freqs[index1] = freqs_line + width * root;
                    nmbrs[index1] = index1;

                    index1++;
                }

                index0++;
            }
        }

        // Sort frequencies
        heapsort (freqs, nmbrs);


        // Set all frequencies nu
        for (Size fl = 0; fl < parameters.nfreqs(); fl++)
        {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }


        // Create lookup table for the frequency corresponding to each line
        Size1 nmbrs_inverted (parameters.nfreqs());

        for (Size fl = 0; fl < parameters.nfreqs(); fl++)
        {
            nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters.nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters.nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters.nfreqs();
        }

        Size index2 = 0;

        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            for (Size k = 0; k < lines.lineProducingSpecies[l].nr_line[p].size(); k++)
            {
                for (Size z = 0; z < lines.lineProducingSpecies[l].nr_line[p][k].size(); z++)
                {
                    lines.lineProducingSpecies[l].nr_line[p][k][z] = nmbrs_inverted[index2];

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
    spectralDiscretisation = SD_Lines;

    return (0);
}


///  Computer for spectral (=frequency) discretisation
///  Gives same frequency bins to each point
///    @param[in] width : corresponding line width for frequency bins
/////////////////////////////////////////////////////////////////////
int Model :: compute_spectral_discretisation (const Real width)
{
    cout << "Computing spectral discretisation..." << endl;

    threaded_for (p, parameters.npoints(),
    {
        Real1 freqs (parameters.nfreqs());
        Size1 nmbrs (parameters.nfreqs());

        Size index0 = 0;//temp index for line number
        Size index1 = 0;//temp index for frequencies


        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real freqs_line = lines.line[index0];

                for (Size z = 0; z < parameters.nquads(); z++)
                {
                    const Real root = lines.lineProducingSpecies[l].quadrature.roots[z];

                    freqs[index1] = freqs_line + freqs_line * width * root;
                    nmbrs[index1] = index1;

                    index1++;
                }

                index0++;
            }
        }


        // Sort frequencies
        heapsort (freqs, nmbrs);


        // Set all frequencies nu
        for (Size fl = 0; fl < parameters.nfreqs(); fl++)
        {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }



        // Create lookup table for the frequency corresponding to each line
        Size1 nmbrs_inverted (parameters.nfreqs());

        for (Size fl = 0; fl < parameters.nfreqs(); fl++)
        {
            nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters.nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters.nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters.nfreqs();
        }

        Size index2 = 0;//temp index for lookup table from above

        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            for (Size k = 0; k < lines.lineProducingSpecies[l].nr_line[p].size(); k++)
            {
                for (Size z = 0; z < lines.lineProducingSpecies[l].nr_line[p][k].size(); z++)
                {
                    lines.lineProducingSpecies[l].nr_line[p][k][z] = nmbrs_inverted[index2];

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


///  Computer for spectral (=frequency) discretisation
///  Gives same frequency bins to each point
///    @param[in] nu_min : minimal frequency
///    @param[in] nu_max : maximal frequency
///////////////////////////////////////////////////////
int Model :: compute_spectral_discretisation (
    const long double nu_min,
    const long double nu_max )
{
    cout << "Computing spectral discretisation..." << endl;

    const long double dnu = (nu_max - nu_min) / (parameters.nfreqs() - 1);

    threaded_for (p, parameters.npoints(),
    {
        for (Size f = 0; f < parameters.nfreqs(); f++)
        {
            radiation.frequencies.nu(p, f) = (Real) (nu_min + f*dnu);

            radiation.frequencies.appears_in_line_integral[f] = false;;
            radiation.frequencies.corresponding_l_for_spec[f] = parameters.nfreqs();
            radiation.frequencies.corresponding_k_for_tran[f] = parameters.nfreqs();
            radiation.frequencies.corresponding_z_for_line[f] = parameters.nfreqs();
        }
    })

    // Set spectral discretisation setting
    spectralDiscretisation = SD_Image;

    return (0);
}


/// Computer for the level populations, assuming LTE
////////////////////////////////////////////////////
int Model :: compute_LTE_level_populations ()
{
    cout << "Computing LTE level populations..." << endl;

    // Initialize levels, emissivities and opacities with LTE values
    lines.iteration_using_LTE (chemistry.species.abundance, thermodynamics.temperature.gas);

    return (0);
}

///  Restarting from the iteration levelpops, assuming it has been written to disk
///    @param[in] iteration : The iteration to start from
///    @param[in] lvl : The coarsening level to start from
//////////////////////////////////////////////////////////////////////////////////
///  Note: when the level populations cannot be read, starts from LTE instead without warning
///  Note: currently no state of where we are in any multilevel operation is stored, so unless changed, please use this only without multiresolution or just naive multiresolution (then the number of iterations done that is reported will be incorrect (because it doesnt count the iterations prior to loading))
int Model :: restart_from_iteration(Size iteration, Size lvl)
{// TODO: currently, the mrController information is NOT SAVED, so this is not that useful for restarting a multiresolution scheme
    compute_LTE_level_populations();
    IoPython io = IoPython ("hdf5", parameters.model_name());
    lines.read_populations_of_iteration(io, iteration, lvl);
    std::cout<<"Read populations from disk"<<std::endl;
    lines.set_emissivity_and_opacity ();
    std::cout<<"Restarting from iteration: "<<iteration<<std::endl;
    iteration_to_start_from=iteration;
    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_shortchar_order_0 ()
{
    cout << "Computing radiation field..." << endl;

    // const Size length_max = 4*parameters.npoints() + 1;
    // const Size  width_max =   parameters.nfreqs ();

    // Solver solver (length_max, width_max, parameters.n_off_diag);

    Solver solver;
    solver.setup <CoMoving>        (*this);
    solver.solve_shortchar_order_0 (*this);

    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_feautrier_order_2 ()
{
    cout << "Computing radiation field..." << endl;

    // const Size length_max = 4*parameters.npoints() + 1;
    // const Size  width_max =   parameters.nfreqs ();

    // Solver solver (length_max, width_max, parameters.n_off_diag);

    Solver solver;

    //ray setup stuff should also be done using the same grid as the actual calculations
    if(using_same_grid)//temporarily do as if we are using the original grid, of course only computing for the lower grid
    {
        geometry.points.multiscale.temporary_set_level_to_original_grid();
    }

    solver.setup <CoMoving>        (*this);

    if(using_same_grid)//temporarily do as if we are using the original grid, of course only computing for the lower grid
    {
        geometry.points.multiscale.reset_temporary_level();
    }


    solver.solve_feautrier_order_2 (*this);

    return (0);
}

// //// also no declaration here in .hpp file; thus commented out
// ///  Computer for the radiation field
// /////////////////////////////////////
// int Model :: compute_radiation_field_2nd_order_Feautrier ()
// {
//     cout << "Computing radiation field..." << endl;
//
//     const Size length_max = 4*parameters.npoints() + 1;
//     const Size  width_max =   parameters.nfreqs ();
//
//
//     cout << "npoints = " << parameters.npoints() << endl;
//     cout << "nfreqs  = " << parameters.nfreqs () << endl;
//
//     cout << "l_max = " << length_max << endl;
//     cout << "w_max = " <<  width_max << endl;
//
//     Solver solver (length_max, width_max, parameters.n_off_diag);
//     solver.solve_2nd_order_Feautrier (*this);
//
//     return (0);
// }


///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model :: compute_Jeff ()
{

    //as usual, only do this for the points currently in the grid
    vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
    Size nbpoints=points_in_grid.size();


    for (LineProducingSpecies &lspec : lines.lineProducingSpecies)
    {
       // Lambda = MatrixXd::Zero (lspec.population.size(), lspec.population.size());

        threaded_for (idx, nbpoints,
        {
            const Size p=points_in_grid[idx];

            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                const Size1 freq_nrs = lspec.nr_line[p][k];

                // Initialize values
                lspec.Jlin[p][k] = 0.0;

                // Integrate over the line
                for (Size z = 0; z < parameters.nquads(); z++)
                {
                    lspec.Jlin[p][k] += lspec.quadrature.weights[z] * radiation.J(p, freq_nrs[z]);
                }


                double diff = 0.0;

                // Collect the approximated part
                for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
                {
                    const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);

                    diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
                }

                lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
                lspec.Jdif[p][k] = HH_OVER_FOUR_PI * diff;
            }
        })
    }

    return (0);
}


///  Compute level populations from statistical equilibrium
///////////////////////////////////////////////////////////
int Model :: compute_level_populations_from_stateq ()
{
    vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();

    lines.iteration_using_statistical_equilibrium (
            chemistry.species.abundance,
            thermodynamics.temperature.gas,
            parameters.pop_prec(),
            points_in_grid                       );

    return (0);
}


///  Compute level populations self-consistenly with the radiation field using multiresolution
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @return Total number of iterations done
///////////////////////////////////////////////////////////////////////////////
/// Assumes one has setup the multiresolution beforehand see 'model.tpp' 'setup_multiresolution'
int Model :: compute_level_populations_multiresolution (
    const bool use_Ng_acceleration)
{
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines)
    {
        throw std::runtime_error ("Spectral discretisation was not set for Lines!");
    }

    //set curr coarsening level to max
    Size max_coars_lvl=geometry.points.multiscale.get_max_coars_lvl();

    int iteration_sum    = 0;//sum of total number of iterations; will be returned when finished

    //storing the number of iterations done at each level
    vector<Size> iterations_per_level;
    iterations_per_level.resize(max_coars_lvl+1);

    //Timing info per multiresolution level
    vector<Timer> timers_per_level;
    for (Size i=0; i<=max_coars_lvl; i++)
    {
        Timer temp_timer("Time on level "+std::to_string(i));
        timers_per_level.push_back(temp_timer);
    }
    Timer totalTime("Total time");
    Timer interpolation_timer("Interpolation time");
    Timer* current_timer_pointer;


    current_timer_pointer=&timers_per_level[max_coars_lvl];
    //For now, I assume we start iterating at the coarsest level and start both timers at the same time
    (*current_timer_pointer).start();
    totalTime.start();

    bool finished=false;//whether the multiresolution procedure has finished
    bool multiresolution_operation_happened=false;//denotes whether we have just done a multiresolution operation: either 'restrict' or 'interpolate_corrections'
    //neccessary because we should calculate convergence compared to previous iteration on SAME GRID LEVEL (when possible).

    // Initialize the number of iterations
    int iteration        = iteration_to_start_from;//when resuming from a saved iteration, this can be greater than 0

    //FIXME: also try to implement a way to also read the previous level pops
    int iteration_normal = 0;// iteration counter for determining when to use Ng-acceleration
    // Also initialize the previous fraction of converged points
    vector<double> prev_it_frac_not_converged(parameters.nlspecs(),1);

    // Initialize errors
    error_mean.clear ();
    error_max .clear ();
    // Initialize some_not_converged; denoting that some level population have not yet converged
    bool some_not_converged = true;

    while(!finished)
    {

        // Iterate as long as some levels are not converged
        switch (mrControllerHelper.get_next_action()) {

        //Signal to go to the coarsest level
        case MrController::Actions::goto_coarsest:
        {//Not much is happening here, so no timer adjustements
            std::cout<<"Action goto_coarsest"<<std::endl;
            geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
            std::cout<<"Current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<std::endl;

            break;
        }

        //Signal to do one iteration the current level
        case MrController::Actions::stay:
        {
            //Stopping current timer and replace with timer at current level
            (*current_timer_pointer).stop();
            timers_per_level[geometry.points.multiscale.get_curr_coars_lvl()].start();
            current_timer_pointer=&timers_per_level[geometry.points.multiscale.get_curr_coars_lvl()];

            std::cout<<"Action stay"<<std::endl;

            iterations_per_level[geometry.points.multiscale.get_curr_coars_lvl()]++;
            iteration++;
            // logger.write ("Starting iteration ", iteration);
            cout << "Starting iteration " << iterations_per_level[geometry.points.multiscale.get_curr_coars_lvl()] << " for coarsening level: " << geometry.points.multiscale.get_curr_coars_lvl() << endl;


            // If a multiresolution operation has happened, just put the previous level populations here to compare against
            if (multiresolution_operation_happened)
            {
                multiresolution_operation_happened=false;
                lines.set_all_level_pops(computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]);
            }

            // Start assuming convergence
            some_not_converged = false;

            if (use_Ng_acceleration && (iteration_normal == 4))
            {
                vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
                //In this next line, we will also calculate the convergence of the levelpops
                lines.iteration_using_Ng_acceleration (parameters.pop_prec(), points_in_grid);

                iteration_normal = 0;
            }
            else
            {
                cout << "Computing the radiation field..." << endl;
                // calculate the radiation field on this level
                compute_radiation_field_feautrier_order_2 ();

                std::cout << "Computed feautrier" << std::endl;
                compute_Jeff                              ();

                vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
                //In this next line, we will also calculate the convergence of the levelpops (directly after computing the level populations)
                lines.iteration_using_statistical_equilibrium (
                    chemistry.species.abundance,
                    thermodynamics.temperature.gas,
                    parameters.pop_prec(),
                    points_in_grid);

                iteration_normal++;
            }


            if (using_same_grid)//if using the same grid, we need to interpolate for the other non-calculated points
            {
                Size current_lvl=mrControllerHelper.get_current_level();
                std::cout<<"current level"<<current_lvl<<std::endl;
                //now interpolate for all other levels
                for (Size level_diff=0; level_diff<current_lvl; level_diff++)
                {
                    interpolate_levelpops_local(current_lvl-level_diff);
                }

                // Again, do not forget to set emissivity and opacity after interpolating the levelpops
                lines.set_emissivity_and_opacity ();

            }


            //If enabled, we now write the level populations to the hdf5 file
            if (parameters.writing_populations_to_disk){
                IoPython io = IoPython ("hdf5", parameters.model_name());
                lines.write_populations_of_iteration(io, iterations_per_level[geometry.points.multiscale.get_curr_coars_lvl()], geometry.points.multiscale.get_curr_coars_lvl());
                std::cout<<"Wrote populations to disk"<<std::endl;
            }

            //And now we check for all species whether the level populations have already converged
            for (int l = 0; l < parameters.nlspecs(); l++)
            {
                //note: this information is only based on the points currently in the grid.
                error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
                error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

                if (lines.lineProducingSpecies[l].fraction_not_converged > 0.005)
                {
                    some_not_converged = true;
                }

                const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;
                prev_it_frac_not_converged[l]=fnc;

                cout << "Already " << 100 * (1.0 - fnc) << " % converged!" << endl;
            }

            //If all level populations have converged, let it know to the mrController
            if (!some_not_converged)
            {
                iteration_sum+=iteration;
                iteration=0;
                //signal convergence on current grid
                mrControllerHelper.converged_on_current_grid();
                std::cout<<"Signaling convergence"<<std::endl;
            }

            break;
        }

        //Signal when one has to interpolate the level populations
        case MrController::Actions::interpolate_levelpops:
        {
            //Stopping current timer and replace with interpolation timer
            (*current_timer_pointer).stop();
            interpolation_timer.start();
            current_timer_pointer=&interpolation_timer;

            std::cout<<"Action interpolate levelpops"<<std::endl;

            //Saving the current level populations
            computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();

            //Now interpolating the level populations to the next finer level
            cout<<"trying to interpolate level populations; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
            interpolate_levelpops_local(geometry.points.multiscale.get_curr_coars_lvl());
            cout<<"successfully interpolated level populations"<<endl;
            geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);

            //after interpolating the level populations, do NOT forget to set the opacities and emissivities, otherwise this was all for nothing...
            lines.set_emissivity_and_opacity ();
            //and also reset some statistics
            // reinitialize the number of iterations
            iteration        = 0;
            iteration_normal = 0;//TODO: figure out if we still use iteration_normal

            // Also reinitialize the previous fraction of converged points
            prev_it_frac_not_converged=vector<double>(parameters.nlspecs(),1);

            // reinitialize errors
            error_mean.clear ();
            error_max .clear ();
            // reinitialize some_not_converged
            some_not_converged = true;

            break;
        }

        //Signals the restriction operator
        //It turns out that, because of the choice of our restriction operator, (mere insertion at points still in grid)
        // that nothing particularly interesting happens to the right-hand side with this restriction operator.
        case MrController::Actions::restrict:
        {
            //Not much is happening here, so no timer stuff
            std::cout<<"Action restrict"<<std::endl;
            multiresolution_operation_happened=true;

            //Saving the current level populations
            computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();

            //Using one step coarser grid
            geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()+1);

            //Normally after restricting the level populations, we should not forget to set the opacities and emissivities
            //However, this is not necessary due to choice of restriction operator
            // lines.set_emissivity_and_opacity ();

            //Also reset some statistics
            // reinitialize the number of iterations
            iteration        = 0;
            iteration_normal = 0;

            // Also initialize the previous fraction of converged points
            prev_it_frac_not_converged=vector<double>(parameters.nlspecs(),1);

            // reinitialize errors
            error_mean.clear ();
            error_max .clear ();
            // reinitialize some_not_converged
            some_not_converged = true;

            break;
        }

        //Signals the interpolation of corrections to the level populations
        //In contrast to the restriction action, this is a bit more interesting. Instead of only interpolating the solution on the coarser grid,
        //  I will also need to restrict and then interpolate an older solution on the finer grid. Then we take the difference of this result and add it to the previous solution (on the finer grid)
        //  For points on the coarser grid, there is will be no difference compared to merely interpolating the coarser level populations;
        //  However, for points in the difference between the two grid, there might be a difference due to the interpolation
        case MrController::Actions::interpolate_corrections:
        {
            //Stopping current timer and replace with interpolation timer
            (*current_timer_pointer).stop();
            interpolation_timer.start();
            current_timer_pointer=&interpolation_timer;

            multiresolution_operation_happened=true;
            std::cout<<"Action interpolate corrections"<<std::endl;

            //In order to make it a bit simpler to implement, we will interpolate the current solution and the restricted solution on the finer grid independently

            //Saving the current level populations
            computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();


            //pseudocode of what will be happening
            //for all species:
            //new levelpops=old levelpops+(interpolated(current levelpops)-interpolated(restricted(old levelpops)))
            vector<VectorXr> old_levelpops=computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()-1];

            //now interpolate current levelpops
            vector<VectorXr> rel_diff_pops;
            rel_diff_pops.resize(parameters.nlspecs());

            for (Size specidx=0; specidx<parameters.nlspecs(); specidx++)
            {//this is just the difference in level populations
                LineProducingSpecies currlspec=lines.lineProducingSpecies[specidx];
                VectorXr currspeclevelpops=computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()][specidx];
                VectorXr diff_pops=currspeclevelpops-old_levelpops[specidx];

                VectorXr temp_rel_diff_pops=diff_pops;//just to get the right size

                //TODO: this way of calculating the corrections is quite slow; try to use operations on the whole vector instead
                for (Size pointidx=0;pointidx<parameters.npoints();pointidx++)
                {
                    if (geometry.points.multiscale.mask[geometry.points.multiscale.get_curr_coars_lvl()][pointidx])//point still in grid, so we can calculate the relative difference
                    {
                        VectorXr temp_inverse_curr_levelpops=currspeclevelpops(Eigen::seq(currlspec.index(pointidx,0),currlspec.index(pointidx,currlspec.linedata.nlev-1))).cwiseInverse();
                        temp_rel_diff_pops(Eigen::seq(currlspec.index(pointidx,0),currlspec.index(pointidx,currlspec.linedata.nlev-1)))=diff_pops(Eigen::seq(currlspec.index(pointidx,0),currlspec.index(pointidx,currlspec.linedata.nlev-1))).cwiseProduct(temp_inverse_curr_levelpops);
                    }
                    else
                    {//otherwise the point does not lie in the coarser grid, thus the relative difference must be equal to zero
                        temp_rel_diff_pops(Eigen::seq(currlspec.index(pointidx,0),currlspec.index(pointidx,currlspec.linedata.nlev-1))).setZero();
                    }
                }

                rel_diff_pops[specidx]=temp_rel_diff_pops;
            }

            //Now finally interpolate the relative differences
            interpolate_relative_differences_local(geometry.points.multiscale.get_curr_coars_lvl(), rel_diff_pops);

            geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);

            vector<VectorXr> corrected_levelpops;
            //TODO: or it is possible that this thing here is slow
            //Now applying the correction and then renormalizing the level populations
            for (Size specidx=0; specidx<parameters.nlspecs(); specidx++)
            {
                // VectorXr temp_corrected_levelpops=old_levelpops[specidx]+interpolated_new_levelpops[specidx]-interpolated_old_levelpops[specidx];
                VectorXr temp_corrected_levelpops=old_levelpops[specidx]+old_levelpops[specidx].cwiseProduct(rel_diff_pops[specidx]);

                //setting negative levelpops to 0.
                temp_corrected_levelpops.cwiseMax(0);

                // For finding out which abundance corresponds to to the current species
                Size speciesnum=lines.lineProducingSpecies[specidx].linedata.num;

                vector<Size> current_points_in_grid=geometry.points.multiscale.get_current_points_in_grid();

                //Note: we might have some minor performance issues with multiple threads accessing data near eachother
                //for every point in the current grid, we will check if some levelpopulations of that point are negative
                //Note: because we are interpolating the relative corrections, we will need to renormalize the levelpops
                threaded_for (p, current_points_in_grid.size(),
                    Size current_point=current_points_in_grid[p];

                    Real abund=static_cast<Real>(chemistry.species.abundance[current_point][speciesnum]);
                    Size segment_start=lines.lineProducingSpecies[specidx].index(current_point,0);
                    Size n_levels=lines.lineProducingSpecies[specidx].linedata.nlev;
                    //Note: because the sum of levelpops must always be equal to the abundance, at least one value should be greater than zero (after zeroing the negative values)
                    // and thus we will never divide by zero
                    Real sum_of_levelpops=temp_corrected_levelpops.segment(segment_start,n_levels).sum();
                    //and finally renormalizing it; such that the sum of levelpops is again equal to the abundance
                    temp_corrected_levelpops.segment(segment_start,n_levels)=temp_corrected_levelpops.segment(segment_start,n_levels)*(abund/sum_of_levelpops);

                )//end of threaded_for
                //and add the corrected levelpops for this species
                corrected_levelpops.push_back(temp_corrected_levelpops);
            }

            std::cout<<"Now setting the corrected levelpops"<<std::endl;
            lines.set_all_level_pops(corrected_levelpops);

            //after interpolating the level populations, do NOT forget to set the opacities and emissivities, otherwise this was all for nothing...
            lines.set_emissivity_and_opacity ();
            //and also reset some statistics
            // reinitialize the number of iterations
            iteration        = 0;
            iteration_normal = 0;
            // Also initialize the previous fraction of converged points
            prev_it_frac_not_converged=vector<double>(parameters.nlspecs(),1);

            // reinitialize errors
            error_mean.clear ();
            error_max .clear ();
            // reinitialize some_not_converged
            some_not_converged = true;

            break;
        }

        //Signals that the multiresolution procedure has finished
        case MrController::Actions::finish:
        {
            //Stop the timers
            (*current_timer_pointer).stop();
            totalTime.stop();

            //And print some statistics
            std::cout<<"Action finish"<<std::endl;
            finished=true;
            std::cout<<"multiresolution sequence is finished"<<std::endl;
            std::cout<<"Nb iterations per level:"<<std::endl;
            for (Size levelidx=0;levelidx<iterations_per_level.size();levelidx++)
            {
                std::cout<<"Level "<<levelidx<<" #iterations: "<<iterations_per_level[levelidx]<<std::endl;
            }

            //Also print out the timings
            std::cout<<totalTime.get_print_total_string()<<std::endl;
            for (Size i=0; i<=max_coars_lvl; i++)
            {
                std::cout<<timers_per_level[i].get_print_total_string()<<std::endl;
            }
            std::cout<<interpolation_timer.get_print_total_string()<<std::endl;

            break;
        }

        //Signals to do nothing, can help with simplifying the multiresolution controller
        case MrController::Actions::do_nothing:
        {
            std::cout<<"Action do nothing"<<std::endl;
            break;
        }

        default:
        {
            std::cout<<"Hey, you shouldn't be here"<<std::endl;
            break;
        }

        }//end of switch statement

        // std::cout<<"here1"<<std::endl;
        // if (using_same_grid)//if using the same grid, we need to interpolate for the other non-calculated points
        // {
        //     Size current_lvl=mrControllerHelper.get_current_level();
        //     std::cout<<"current level"<<current_lvl<<std::endl;
        //     //now interpolate for all other levels
        //     for (Size level_diff=0; level_diff<current_lvl; level_diff++)
        //     {
        //         interpolate_levelpops_local(current_lvl-level_diff);
        //     }
        // }

    }//end of giant while loop until finished

    return iteration_sum;
}


///  Compute level populations self-consistenly with the radiation field
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return Number of iterations done
///////////////////////////////////////////////////////////////////////////////
int Model :: compute_level_populations (
    const bool use_Ng_acceleration,
    const long max_niterations     )
{
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines)
    {
        throw std::runtime_error ("Spectral discretisation was not set for Lines!");
    }

    // Initialize the number of iterations
    int iteration        = iteration_to_start_from;
    int iteration_normal = 0;

    // Initialize errors
    error_mean.clear ();
    error_max .clear ();

    // Initialize some_not_converged
    bool some_not_converged = true;

    // Iterate as long as some levels are not converged
    while (some_not_converged && (iteration < max_niterations))
    {
        iteration++;

        cout << "Starting iteration " << iteration << endl;

        // Start assuming convergence
        some_not_converged = false;

        if (use_Ng_acceleration && (iteration_normal == 4))
        {
            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
            lines.iteration_using_Ng_acceleration (parameters.pop_prec(), points_in_grid);

            iteration_normal = 0;
        }
        else
        {
            cout << "Computing the radiation field..." << endl;

            compute_radiation_field_feautrier_order_2 ();
            compute_Jeff                              ();

            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();

            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec(),
                points_in_grid);

            iteration_normal++;
        }

        //If enabled, we now write the level populations to the hdf5 file
        if (parameters.writing_populations_to_disk){
            IoPython io = IoPython ("hdf5", parameters.model_name());
            lines.write_populations_of_iteration(io, iteration, geometry.points.multiscale.get_curr_coars_lvl());
            std::cout<<"Wrote populations to disk"<<std::endl;
        }


        for (int l = 0; l < parameters.nlspecs(); l++)
        {
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            if (lines.lineProducingSpecies[l].fraction_not_converged > 0.005)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

            cout << "Already " << 100 * (1.0 - fnc) << " % converged!" << endl;
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;

    return iteration;
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_image (const Size ray_nr)
{
    cout << "Computing image..." << endl;

    // const Size length_max = 4*parameters.npoints() + 1;
    // const Size  width_max =   parameters.nfreqs ();

    // Solver solver (length_max, width_max, parameters.n_off_diag);

    Solver solver;
    solver.setup <Rest>            (*this);
    solver.image_feautrier_order_2 (*this, ray_nr);

    return (0);
}

/// Sets the emmisivity and opacity
///////////////////////////////////
int Model :: set_eta_and_chi ()
{
    Solver solver;
    solver.set_eta_and_chi (*this);

    return (0);
}

/// Sets the boundary conditions
////////////////////////////////
int Model :: set_boundary_condition ()
{
    Solver solver;
    solver.set_boundary_condition (*this);

    return (0);
}
