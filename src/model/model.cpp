#include "paracabs.hpp"
#include "model.hpp"
#include "tools/heapsort.hpp"
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
    // geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
    int iteration_sum    = 0;
    // Now testing without the finest grid
    for(Size subtract=0; subtract<=max_coars_lvl-1; subtract++)
    {
      Size curr_max_coars_lvl=max_coars_lvl-subtract;
      std::cout<<"Current coarsening level: "<< curr_max_coars_lvl<<std::endl;
      std::cout<<"Current number points:"<<geometry.points.multiscale.get_total_points(curr_max_coars_lvl)<<std::endl;
      geometry.points.multiscale.set_curr_coars_lvl(curr_max_coars_lvl);

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
            lines.iteration_using_Ng_acceleration (parameters.pop_prec());

            iteration_normal = 0;
        }
        else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;
            // calculate radiation field on coarser level
            geometry.points.multiscale.set_curr_coars_lvl(curr_max_coars_lvl);
            compute_radiation_field_feautrier_order_2 ();
            // and interpolate it
            // while(geometry.points.multiscale.get_curr_coars_lvl()>0)
            // {//maybe TODO: add support for interpolating skipping levels
            //
            //   const vector<bool> curr_mask=geometry.points.multiscale.get_mask(geometry.points.multiscale.get_curr_coars_lvl());
            //   for (Size test=0;test<parameters.npoints();test++)
            //   {
            //     if (curr_mask[test])
            //     {
            //     std::cout<<radiation.J(test,0)<<std::endl;
            //     }
            //   }
            //   cout<<"trying to interpolate matrix; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
            //   //TODO print number of nans
            //   //for all frequencies, interpolate J
            //   interpolate_matrix_local(geometry.points.multiscale.get_curr_coars_lvl(),radiation.J);
            //   cout<<"successfully interpolated matrix"<<endl;
            //   geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);
            //   for (Size test=0;test<parameters.npoints();test++)
            //   {
            //     std::cout<<radiation.J(test,0)<<std::endl;
            //   }
            // }
            std::cout << "Computed feautrier" << std::endl;
            compute_Jeff                              ();
            //TODO: also interpolate lambda diagonal thing
            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();
            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec(),
                points_in_grid);

            //Now interpolating the level populations to the finest level//TODO: thus should only happen when we change grid, but does not incurr too high a cost
            while(geometry.points.multiscale.get_curr_coars_lvl()>0)
            {//maybe TODO: add support for interpolating skipping levels
                cout<<"trying to interpolate level populations; current coarsening level: "<<geometry.points.multiscale.get_curr_coars_lvl()<<endl;
                //TODO print number of nans
                //for all frequencies, interpolate J
                interpolate_levelpops_local(geometry.points.multiscale.get_curr_coars_lvl());
                cout<<"successfully interpolated level populations"<<endl;
                geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);
            }

            iteration_normal++;
        }


        for (int l = 0; l < parameters.nlspecs(); l++)
        {
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            // fraction allowed to not be converged:
            const double max_frac_not_converged=(parameters.npoints()-geometry.points.multiscale.get_total_points(curr_max_coars_lvl))/(double)parameters.npoints()+0.005;

            cout << "max frac non coverged: " <<max_frac_not_converged<<endl;

            //check whether the fraction non-corverged points has truly stabilized
            if (((lines.lineProducingSpecies[l].fraction_not_converged > max_frac_not_converged)
                ||(abs(lines.lineProducingSpecies[l].fraction_not_converged-prev_it_frac_not_converged[l])>0.005))
              &&(lines.lineProducingSpecies[l].fraction_not_converged > 0.005))
            // if (lines.lineProducingSpecies[l].fraction_not_converged > 0.005)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;
            prev_it_frac_not_converged[l]=fnc;

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100 * (1.0 - fnc) << " % converged!" << endl;
        }
      } // end of while loop of iterations
      // //TODO add check of final iteration
      // cout<<"trying to interpolate matrix"<<endl;
      // for (Size test=0;test<parameters.npoints();test++)
      // {
      //   std::cout<<radiation.J(test,0)<<std::endl;
      // }
      // //for all frequencies, interpolate J
      // interpolate_matrix_local(geometry.points.multiscale.get_curr_coars_lvl(),radiation.J);
      // cout<<"successfully interpolated matrix"<<endl;
      // geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);

      // Print convergence stats
      cout << "Converged after " << iteration << " iterations" << endl;
      // curr_max_coars_lvl-=1;
      iteration_sum+=iteration;

    }

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
            lines.iteration_using_Ng_acceleration (parameters.pop_prec());

            iteration_normal = 0;
        }
        else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            compute_radiation_field_feautrier_order_2 ();
            cout << "Computed feautrier" << std::endl;
            compute_Jeff                              ();

            vector<Size> points_in_grid=geometry.points.multiscale.get_current_points_in_grid();

            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec(),
                points_in_grid);

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
