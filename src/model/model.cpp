#include "paracabs.hpp"
#include "model.hpp"
#include "tools/heapsort.hpp"
#include "tools/timer.hpp"
#include "solver/solver.hpp"
// #include "mgController/mgControllerHelper.hpp"


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

        Size index0 = 0;
        Size index1 = 0;


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
    spectralDiscretisation = SD_Image;

    return (0);
}


///  Computer for spectral (=frequency) discretisation
///  Gives same frequency bins to each point
///    @param[in] min : minimal frequency
///    @param[in] max : maximal frequency
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
    solver.setup <CoMoving>        (*this);
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


///  compute level populations from statistical equilibrium
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


///TODO remove maxn_iterations; has no place here anymore
///  Compute level populations self-consistenly with the radiation field using multigrid
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] io                  : io object (for writing level populations)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return number of iteration done
///////////////////////////////////////////////////////////////////////////////
int Model :: compute_level_populations_multigrid (
    const bool use_Ng_acceleration,
    const long max_niterations     )
{
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines)
    {
        throw std::runtime_error ("Spectral discretisation was not set for Lines!");
    }

    //set curr coarsening level to max
    Size max_coars_lvl=geometry.points.multiscale.get_max_coars_lvl();
    // std::cout<<"max_coars_lvl: "<<max_coars_lvl<<std::endl;
    // geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
    int iteration_sum    = 0;

    //storing the number of ieterations each level
    vector<Size> iterations_per_level;
    iterations_per_level.resize(max_coars_lvl+1);
    //Timing info
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
    //For now, i assume we start iterating at the coarsest level and start both timers at the same time
    (*current_timer_pointer).start();
    totalTime.start();



    bool finished=false;


    bool multigrid_operation_happened=false;//denotes whether we have just done a multigrid operation: either 'restrict' or 'interpolate_corrections'
    //neccessary because we should calculate convergence compared to previous iteration on SAME GRID LEVEL (when possible).

    //TODO: initialize multigrid

    //TODO: use something a bit more robust and RESET THE CONTROLLER
    // Size curr_max_coars_lvl=mgControllerHelper.get_current_level();//max_coars_lvl-subtract;
    // std::cout<<"Current coarsening level: "<< curr_max_coars_lvl<<std::endl;
    // std::cout<<"Current number points:"<<geometry.points.multiscale.get_total_points(curr_max_coars_lvl)<<std::endl;
    // geometry.points.multiscale.set_curr_coars_lvl(curr_max_coars_lvl);

    // Initialize the number of iterations
    int iteration        = 0;
    int iteration_normal = 0;
    // Also initialize the previous fraction of converged points
    vector<double> prev_it_frac_not_converged(parameters.nlspecs(),1);

    // Initialize errors
    error_mean.clear ();
    error_max .clear ();
    // Initialize some_not_converged
    bool some_not_converged = true;

    while(!finished)
    {
      std::cout<<"looping"<<std::endl;
    // Now testing without the finest grid
    // for(Size subtract=0; subtract<=max_coars_lvl; subtract++)
    // {
      //TODO: rename

      // Iterate as long as some levels are not converged
      switch (mgControllerHelper.get_next_action()) {

        case::MgController::Actions::goto_coarsest:
        {//Not much is happening here, so no timer
        std::cout<<"Action goto_coarsest"<<std::endl;
        geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
        std::cout<<"Current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<std::endl;
        // std::cout<<"Current coarsening level: "<<mgControllerHelper.get_current_level()<<std::endl;

      break;
        }


        case MgController::Actions::stay:
        {
          //Stopping current timer and replace with timer at current level
          (*current_timer_pointer).stop();
          timers_per_level[geometry.points.multiscale.get_curr_coars_lvl()].start();
          current_timer_pointer=&timers_per_level[geometry.points.multiscale.get_curr_coars_lvl()];

        std::cout<<"Action stay"<<std::endl;

      // while (some_not_converged && (iteration < max_niterations))
      // {
        iterations_per_level[geometry.points.multiscale.get_curr_coars_lvl()]++;
        iteration++;
        // logger.write ("Starting iteration ", iteration);
        cout << "Starting iteration " << iterations_per_level[geometry.points.multiscale.get_curr_coars_lvl()] << " for coarsening level: " << geometry.points.multiscale.get_curr_coars_lvl() << endl;


        // If a multigrid operation has happened, just put the previous level populations here to compare against
        if (multigrid_operation_happened)
        {
          multigrid_operation_happened=false;
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
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;
            // calculate radiation field on coarser level
            // geometry.points.multiscale.set_curr_coars_lvl(curr_max_coars_lvl);
            compute_radiation_field_feautrier_order_2 ();

            std::cout << "Computed feautrier" << std::endl;
            compute_Jeff                              ();
            //TODO: also interpolate lambda diagonal thing
            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
            //In this next line, we will also calculate the convergence of the levelpops
            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec(),
                points_in_grid);

            iteration_normal++;
        }


        for (int l = 0; l < parameters.nlspecs(); l++)
        {
            //note: this information is only based on the points currently in the grid.
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            // fraction allowed to not be converged:
            // const double max_frac_not_converged=(parameters.npoints()-geometry.points.multiscale.get_total_points(curr_max_coars_lvl))/(double)parameters.npoints()+0.005;

            // cout << "max frac non coverged: " <<max_frac_not_converged<<endl;

            //check whether the fraction non-corverged points has truly stabilized
            // if (((lines.lineProducingSpecies[l].fraction_not_converged > max_frac_not_converged)
            //     ||(abs(lines.lineProducingSpecies[l].fraction_not_converged-prev_it_frac_not_converged[l])>0.005))
            //   &&(lines.lineProducingSpecies[l].fraction_not_converged > 0.005))
            if (lines.lineProducingSpecies[l].fraction_not_converged > 0.005)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;
            prev_it_frac_not_converged[l]=fnc;

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100 * (1.0 - fnc) << " % converged!" << endl;
        }
      // } // end of while loop of iterations
      // //TODO add check of final iteration
      // std::cout<<"iteration: "<<iteration<<std::endl;
      // std::cout<<"max iteration: "<<max_niterations<<std::endl;
      //if converged //or reached max number of iterations; moved logic to mgController
      if (!some_not_converged)//||iteration==max_niterations)
      {
      // Print convergence stats
      // cout << "Converged after " << iteration << " iterations" << endl;
      // curr_max_coars_lvl-=1;
      iteration_sum+=iteration;
      iteration=0;
      //signal convergence on current grid
      mgControllerHelper.converged_on_current_grid();
      std::cout<<"Signaling convergence"<<std::endl;

      }
      break;
        }

        case MgController::Actions::interpolate_levelpops:
        {
          //Stopping current timer and replace with timer at current level
          (*current_timer_pointer).stop();
          interpolation_timer.start();
          current_timer_pointer=&interpolation_timer;

        std::cout<<"Action interpolate levelpops"<<std::endl;

        //Saving the current level populations
        computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();

      //Now interpolating the level populations to the next finer level
      // if (geometry.points.multiscale.get_curr_coars_lvl()>0)
      {//maybe TODO: add support for interpolating skipping levels
          cout<<"trying to interpolate level populations; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
          //TODO print number of nans
          //for all frequencies, interpolate J
          interpolate_levelpops_local(geometry.points.multiscale.get_curr_coars_lvl());
          cout<<"successfully interpolated level populations"<<endl;
          geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);
      }
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

        //It turns out that, because of the choice of our restriction operator, (mere insertation at points still in grid)
        // that nothing particularly interesting happens to the right-hand side with this restriction operator.
        case MgController::Actions::restrict:
        {
          //Not much is happening here, so no timer

          std::cout<<"Action restrict"<<std::endl;
          multigrid_operation_happened=true;

          //Saving the current level populations
          computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();


          //Using one step coarser grid
          geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()+1);

          //after restricting the level populations, do NOT forget to set the opacities and emissivities, otherwise this was all for nothing...
          // Not necessary to set emissivities and opacities due to choice of restriction Operator
          // lines.set_emissivity_and_opacity ();
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

        //In contrast to the restriction action, this is a bit more interesting. Instead of only interpolating the solution on the coarser grid,
        // i will also need to restrict and then interpolate an older solution on the finer grid. Then we take the difference of this result and add it to the previous solution (on the finer grid)
        // For points on the coarser grid, there is will be no difference compared to merely interpolating the coarser level populations;
        // However, for points in the difference between the two grid, there might be a difference due to the interpolation
        case::MgController::Actions::interpolate_corrections:
        {
          //Stopping current timer and replace with timer at current level
          (*current_timer_pointer).stop();
          interpolation_timer.start();
          current_timer_pointer=&interpolation_timer;


          multigrid_operation_happened=true;
          std::cout<<"Action interpolate corrections"<<std::endl;

          //In order to make it a bit simpler to implement, we will interpolate the current solution and the restricted solution on the finer grid independently

          //Saving the current level populations
          computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()]=lines.get_all_level_pops();

          //for all species:
          //new levelpops=old levelpops+(interpolated(current levelpops)+interpolated(restricted(old levelpops)))
          vector<VectorXr> old_levelpops=computed_level_populations[geometry.points.multiscale.get_curr_coars_lvl()-1];
          //now interpolate current levelpops


          //Maybe TODO: create function that directly interpolates both, such that we do not waste extra time when recreating the same matrices
          // cout<<"trying to interpolate level populations; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
          interpolate_levelpops_local(geometry.points.multiscale.get_curr_coars_lvl());
          // cout<<"successfully interpolated level populations"<<endl;

          vector<VectorXr> interpolated_new_levelpops=lines.get_all_level_pops();

          // cout<<"trying to interpolate old level populations; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
          interpolate_levelpops_local(geometry.points.multiscale.get_curr_coars_lvl());
          // cout<<"successfully interpolated level populations"<<endl;

          vector<VectorXr> interpolated_old_levelpops=lines.get_all_level_pops();

          geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);

          vector<VectorXr> corrected_levelpops;
          //FIXME: add some check for negative values maybe and do something about them
          for (Size specidx=0; specidx<parameters.nlspecs(); specidx++)
          {
            VectorXr temp_corrected_levelpops=old_levelpops[specidx]+interpolated_new_levelpops[specidx]-interpolated_old_levelpops[specidx];

            //setting negative levelpops to 0.
            temp_corrected_levelpops.cwiseMax(0);

            // For finding out which abundance corresponds to to the current species
            Size speciesnum=lines.lineProducingSpecies[specidx].linedata.num;
            // We will need to renormalize the level pops, so first we should collect them
            // vector<Real> linefracs;
            // linefracs.resize(lines.lineProducingSpecies[specidx].linedata.nlev);


            vector<Size> current_points_in_grid=geometry.points.multiscale.get_current_points_in_grid();


            //Note: we might have some issues with multiple threads accessing data near eachother
            //for every point in the current grid, we will check if some levelpopulations of that point are negative
            threaded_for (p, current_points_in_grid.size(),
              Size current_point=current_points_in_grid[p];
              bool negative_value=false;
              for (Size levidx=0; levidx<lines.lineProducingSpecies[specidx].linedata.nlev; levidx++)
              {
                if (temp_corrected_levelpops(lines.lineProducingSpecies[specidx].index(current_point,levidx))==0)
                {
                  //probably, there was some negative value there
                  negative_value=true;
                  break;
                  // //and set the negative value to 0
                  // temp_corrected_levelpops(lines.lineProducingSpecies[specidx].index(current_point,levidx))=0;
                }
              }

              //if negative level population encountered, renormalize
              if (negative_value)
              {
                Real abund=static_cast<Real>(chemistry.species.abundance[current_point][speciesnum]);
                Size segment_start=lines.lineProducingSpecies[specidx].index(current_point,0);
                Size nb_levels=lines.lineProducingSpecies[specidx].linedata.nlev;
                //Note: because the sum of levelpops must always be equal to the abundance, at least one value should be greater than zero (after zeroing the negative values)
                // and thus we will never divide by zero
                Real sum_of_levelpops=temp_corrected_levelpops.segment(segment_start,nb_levels).sum();
                //and finally renormalizing it; such that the sum of levelpops is again equal to the abundance
                temp_corrected_levelpops.segment(segment_start,nb_levels)=temp_corrected_levelpops.segment(segment_start,nb_levels)*(abund/sum_of_levelpops);
              }
            )

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


        case MgController::Actions::finish:
        {
        //Stop the timers
        (*current_timer_pointer).stop();
        totalTime.stop();

        std::cout<<"Action finish"<<std::endl;
        finished=true;
        std::cout<<"Multigrid sequence is finished"<<std::endl;
        std::cout<<"Nb iterations per level:"<<std::endl;
        for (Size levelidx=0;levelidx<iterations_per_level.size();levelidx++)
        {
          std::cout<<"Level "<<levelidx<<" #iterations: "<<iterations_per_level[levelidx]<<std::endl;
        }

        //and print out the timings
        std::cout<<totalTime.get_print_total_string()<<std::endl;
        for (Size i=0; i<=max_coars_lvl; i++)
        {
          std::cout<<timers_per_level[i].get_print_total_string()<<std::endl;
        }
        std::cout<<interpolation_timer.get_print_total_string()<<std::endl;

        }

      break;

        case MgController::Actions::do_nothing:
        {
        std::cout<<"Action do nothing"<<std::endl;
        }

      break;

        default:
        {
        std::cout<<"Hey, you shouldn't be here"<<std::endl;
        break;
        }
      }//end of switch statement
    }//end of giant while loop until finished

    return iteration_sum;
}

// int Model :: compute_level_populations_multigrid (
//     const bool use_Ng_acceleration,
//     const long max_niterations     )
// {
//   //set curr coarsening level to max
//   geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
//   std::cout<<"Current coarsening level: "<< geometry.points.multiscale.get_curr_coars_lvl()<<std::endl;
//   std::cout<<"Current number points:"<<geometry.points.multiscale.get_total_points(geometry.points.multiscale.get_curr_coars_lvl())<<std::endl;
//   //solve and interpolate J for all coarser grids
//   while(geometry.points.multiscale.get_curr_coars_lvl()>0)
//   {
//     compute_level_populations(use_Ng_acceleration,max_niterations);
//     cout<<"trying to interpolate matrix"<<endl;
//     //for all frequencies, interpolate J
//     interpolate_matrix_local(geometry.points.multiscale.get_curr_coars_lvl(),radiation.J);
//     cout<<"successfully interpolated matrix"<<endl;
//     geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);
//     std::cout<<"Current coarsening level: "<< geometry.points.multiscale.get_curr_coars_lvl()<<std::endl;
//     std::cout<<"Current number points:"<<geometry.points.multiscale.get_total_points(geometry.points.multiscale.get_curr_coars_lvl())<<std::endl;
//     // compute_Jeff                              ();
//     // lines.iteration_using_statistical_equilibrium (
//     //     chemistry.species.abundance,
//     //     thermodynamics.temperature.gas,
//     //     parameters.pop_prec()                     );
//   }
//     //finally, solve for the final grid
//     //as a test, we can just not calculate the radiation field again for the finest grid
//     // compute_Jeff                              ();
//     // lines.iteration_using_statistical_equilibrium (
//     //     chemistry.species.abundance,
//     //     thermodynamics.temperature.gas,
//     //     parameters.pop_prec()                     );
//     compute_level_populations(use_Ng_acceleration,max_niterations);
//     return (0);
// }


///  Compute level populations self-consistenly with the radiation field
///  assuming statistical equilibrium (detailed balance for the levels)
///  @param[in] io                  : io object (for writing level populations)
///  @param[in] use_Ng_acceleration : true if Ng acceleration has to be used
///  @param[in] max_niterations     : maximum number of iterations
///  @return number of iteration done
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
    int iteration        = 0;
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

        // logger.write ("Starting iteration ", iteration);
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
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            compute_radiation_field_feautrier_order_2 ();
            cout << "Computed feautrier" << std::endl;
            compute_Jeff                              ();
            cout << "Computed Jeff" << std::endl;

            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();

            cout<< "got point_in_grid"<<std::endl;
            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec(),
                points_in_grid);

              cout<<"iterated using statistical equilibrium"<<std::endl;

            iteration_normal++;
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

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100 * (1.0 - fnc) << " % converged!" << endl;
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;
    cout << "Some extra output" <<endl;

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


int Model :: set_eta_and_chi ()
{
    Solver solver;
    solver.set_eta_and_chi (*this);

    return (0);
}


int Model :: set_boundary_condition ()
{
    Solver solver;
    solver.set_boundary_condition (*this);

    return (0);
}
