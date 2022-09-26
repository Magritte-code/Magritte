#include "paracabs.hpp"
#include "model.hpp"
#include "tools/heapsort.hpp"
#include "solver/solver.hpp"
#include "cooling/cooling.hpp"


void Model :: read (const Io& io)
{
    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Reading Model...                         " << endl;
    cout << "-------------------------------------------" << endl;
    cout << " model file = " << io.io_file                << endl;
    cout << "-------------------------------------------" << endl;

    parameters   ->read (io);
    geometry      .read (io);
    chemistry     .read (io);
    thermodynamics.read (io);
    lines         .read (io);
    radiation     .read (io);
    cooling       .read (io);

    cout << "                                           " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  Model read, parameters:                  " << endl;
    cout << "-------------------------------------------" << endl;
    cout << "  npoints    = " << parameters->npoints   () << endl;
    cout << "  nrays      = " << parameters->nrays     () << endl;
    cout << "  nboundary  = " << parameters->nboundary () << endl;
    cout << "  nfreqs     = " << parameters->nfreqs    () << endl;
    cout << "  nspecs     = " << parameters->nspecs    () << endl;
    cout << "  nlspecs    = " << parameters->nlspecs   () << endl;
    cout << "  nlines     = " << parameters->nlines    () << endl;
    cout << "  nquads     = " << parameters->nquads    () << endl;
    cout << "-------------------------------------------" << endl;
    cout << "                                           " << endl;
}


void Model :: write (const Io& io) const
{
    // Let only root (rank 0) process write output
    if (pc::message_passing::comm_rank() == 0)
    {
        parameters   ->write (io);
        geometry      .write (io);
        chemistry     .write (io);
        thermodynamics.write (io);
        lines         .write (io);
        radiation     .write (io);
        cooling       .write (io);
    }
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

    threaded_for (p, parameters->npoints(),
    {
        Real1 freqs (parameters->nfreqs());
        Size1 nmbrs (parameters->nfreqs());

        Size index0 = 0;
        Size index1 = 0;

        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            const Real inverse_mass = lines.lineProducingSpecies[l].linedata.inverse_mass;

            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real freqs_line = lines.line[index0];
                const Real width      = freqs_line * thermodynamics.profile_width (inverse_mass, p);

                for (Size z = 0; z < parameters->nquads(); z++)
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
        for (Size fl = 0; fl < parameters->nfreqs(); fl++)
        {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }


        // Create lookup table for the frequency corresponding to each line
        Size1 nmbrs_inverted (parameters->nfreqs());

        for (Size fl = 0; fl < parameters->nfreqs(); fl++)
        {
            nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters->nfreqs();
        }

        Size index2 = 0;

        for (Size l = 0; l < parameters->nlspecs(); l++)
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


    // Set corresponding frequencies
    for (Size f = 0; f < parameters->nfreqs(); f++)
    {
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
int Model :: compute_spectral_discretisation (const Real width)
{
    cout << "Computing spectral discretisation..." << endl;

    threaded_for (p, parameters->npoints(),
    {
        Real1 freqs (parameters->nfreqs());
        Size1 nmbrs (parameters->nfreqs());

        Size index0 = 0;
        Size index1 = 0;


        // Add the line frequencies (over the profile)
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lines.lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real freqs_line = lines.line[index0];

                for (Size z = 0; z < parameters->nquads(); z++)
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
        for (Size fl = 0; fl < parameters->nfreqs(); fl++)
        {
            radiation.frequencies.nu(p, fl) = freqs[fl];
        }



        // Create lookup table for the frequency corresponding to each line
        Size1 nmbrs_inverted (parameters->nfreqs());

        for (Size fl = 0; fl < parameters->nfreqs(); fl++)
        {
            nmbrs_inverted[nmbrs[fl]] = fl;

            radiation.frequencies.appears_in_line_integral[fl] = false;;
            radiation.frequencies.corresponding_l_for_spec[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[fl] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[fl] = parameters->nfreqs();
        }

        Size index2 = 0;

        for (Size l = 0; l < parameters->nlspecs(); l++)
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

    const long double dnu = (nu_max - nu_min) / (parameters->nfreqs() - 1);

    threaded_for (p, parameters->npoints(),
    {
        for (Size f = 0; f < parameters->nfreqs(); f++)
        {
            radiation.frequencies.nu(p, f) = (Real) (nu_min + f*dnu);

            radiation.frequencies.appears_in_line_integral[f] = false;;
            radiation.frequencies.corresponding_l_for_spec[f] = parameters->nfreqs();
            radiation.frequencies.corresponding_k_for_tran[f] = parameters->nfreqs();
            radiation.frequencies.corresponding_z_for_line[f] = parameters->nfreqs();
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

    Solver solver;
    solver.setup <CoMoving>        (*this);
    solver.solve_shortchar_order_0 (*this);

    return (0);
}


/// NOTE TO SELF: I might need to template this up if I also add radiative pressure
///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_feautrier_order_2 ()
{
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup <CoMoving> (*this);

    if (parameters->one_line_approximation)
    {
        solver.solve_feautrier_order_2 <OneLine> (*this);
    }
    else
    {
        solver.solve_feautrier_order_2 <None> (*this);
    }

    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_feautrier_order_2_uv ()
{
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup <CoMoving>                  (*this);

    if (parameters->one_line_approximation)
    {
        solver.solve_feautrier_order_2_uv <OneLine> (*this);
    }
    else
    {
        solver.solve_feautrier_order_2_uv <None> (*this);
    }

    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_feautrier_order_2_anis ()
{
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup <CoMoving>                    (*this);

    if (parameters->one_line_approximation)
    {
        solver.solve_feautrier_order_2_anis <OneLine> (*this);
    }
    else
    {
        solver.solve_feautrier_order_2_anis <None> (*this);
    }

    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_feautrier_order_2_sparse ()
{
    cout << "Computing radiation field..." << endl;

    Solver solver;
    solver.setup <CoMoving>                      (*this);

    if (parameters->one_line_approximation)
    {
        solver.solve_feautrier_order_2_sparse <OneLine> (*this);
    }
    else
    {
        solver.solve_feautrier_order_2_sparse <None> (*this);
    }

    return (0);
}


///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model :: compute_Jeff ()
{
    for (LineProducingSpecies &lspec : lines.lineProducingSpecies)
    {
        threaded_for (p, parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                const Size1 freq_nrs = lspec.nr_line[p][k];

                // Initialize values
                lspec.Jlin[p][k] = 0.0;

                // Integrate over the line
                for (Size z = 0; z < parameters->nquads(); z++)
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


///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model :: compute_Jeff_sparse ()
{
    for (LineProducingSpecies &lspec : lines.lineProducingSpecies)
    {
        threaded_for (p, parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                double diff = 0.0;

                // Collect the approximated part
                for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
                {
                    const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);

                    diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
                }

                lspec.Jlin[p][k] = lspec.J(p,k);
                lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
            }
        })
    }

    return (0);
}


///  compute level populations from statistical equilibrium
///////////////////////////////////////////////////////////
int Model :: compute_level_populations_from_stateq ()
{
    lines.iteration_using_statistical_equilibrium (
            chemistry.species.abundance,
            thermodynamics.temperature.gas,
            parameters->pop_prec                 );

    return (0);
}


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
            lines.iteration_using_Ng_acceleration (parameters->pop_prec);

            iteration_normal = 0;
        }
        else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_feautrier_order_2 ();
            compute_Jeff                              ();

            timer_1.stop();
            timer_1.print_total();


            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters->pop_prec                     );

            timer_2.stop();
            timer_2.print_total();


            iteration_normal++;
        }


        for (int l = 0; l < parameters->nlspecs(); l++)
        {
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            if (lines.lineProducingSpecies[l].fraction_not_converged > 1.0 - parameters->convergence_fraction)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
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
int Model :: compute_level_populations_sparse (
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
            lines.iteration_using_Ng_acceleration (parameters->pop_prec);

            iteration_normal = 0;
        }
        else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_feautrier_order_2_sparse ();
            compute_Jeff_sparse                              ();

            timer_1.stop();
            timer_1.print_total();


            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium_sparse (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters->pop_prec                            );

            timer_2.stop();
            timer_2.print_total();


            iteration_normal++;
        }


        for (int l = 0; l < parameters->nlspecs(); l++)
        {
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            if (lines.lineProducingSpecies[l].fraction_not_converged > 1.0 - parameters->convergence_fraction)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
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
int Model :: compute_level_populations_shortchar (
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
            lines.iteration_using_Ng_acceleration (parameters->pop_prec);

            iteration_normal = 0;
        }
        else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            Timer timer_1("Compute Radiation Field");
            timer_1.start();

            compute_radiation_field_shortchar_order_0 ();
            compute_Jeff                              ();

            timer_1.stop();
            timer_1.print_total();


            Timer timer_2("Compute Statistical Equilibrium");
            timer_2.start();

            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters->pop_prec                            );

            timer_2.stop();
            timer_2.print_total();


            iteration_normal++;
        }


        for (int l = 0; l < parameters->nlspecs(); l++)
        {
            error_mean.push_back (lines.lineProducingSpecies[l].relative_change_mean);
            error_max .push_back (lines.lineProducingSpecies[l].relative_change_max);

            if (lines.lineProducingSpecies[l].fraction_not_converged > 1.0 - parameters->convergence_fraction)
            {
                some_not_converged = true;
            }

            const double fnc = lines.lineProducingSpecies[l].fraction_not_converged;

            // logger.write ("Already ", 100 * (1.0 - fnc), " % converged!");
            cout << "Already " << 100.0 * (1.0 - fnc) << " % converged!" << endl;
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

    Solver solver;
    solver.setup <Rest>            (*this);
    solver.image_feautrier_order_2 (*this, ray_nr);

    return (0);
}


///  Computer for image in one point
////////////////////////////////////
int Model :: compute_image_for_point (const Size ray_nr, const Size p)
{
    cout << "Computing image for point " << p << "..." << endl;

    Solver solver;
    solver.setup <Rest>                      (*this);
    solver.image_feautrier_order_2_for_point (*this, ray_nr, p);

    return (0);
}


///  Computer for optical depth iimage
//////////////////////////////////////
int Model :: compute_image_optical_depth (const Size ray_nr)
{
    Solver solver;
    solver.setup <Rest>        (*this);
    solver.image_optical_depth (*this, ray_nr);

    return (0);
}

int Model :: compute_cooling_collisional ()
{
    Vector<Real>& gas_temperature=thermodynamics.temperature.gas;
    Double2 abundance=chemistry.species.abundance;

    lines.compute_cooling_collisional(cooling, abundance, gas_temperature);

    return (0);
}

int Model :: compute_cooling_radiative ()
{
      Solver solver;
      solver.setup <CoMoving> (*this);
      const bool IS_SPARSE=true;
      const bool COMPUTE_UV=false;
      const bool COMPUTE_ANIS=false;
      const bool COMPUTE_LAMBDA=false;
      const bool COMPUTE_COOLING=true;

      if (parameters->one_line_approximation)
      {
          solver.solve_feautrier_order_2 <OneLine, IS_SPARSE, COMPUTE_UV, COMPUTE_ANIS, COMPUTE_LAMBDA, COMPUTE_COOLING> (*this);
      }
      else
      {
          solver.solve_feautrier_order_2 <None, IS_SPARSE, COMPUTE_UV, COMPUTE_ANIS, COMPUTE_LAMBDA, COMPUTE_COOLING> (*this);
      }

      return (0);
}


int Model :: set_eta_and_chi (const Size rr)
{
    Solver solver;
    solver.set_eta_and_chi (*this, rr);

    return (0);
}


int Model :: set_boundary_condition ()
{
    Solver solver;
    solver.set_boundary_condition (*this);

    return (0);
}


int Model :: set_column ()
{
    Solver solver;
    solver.set_column (*this);

    return (0);
}


// /  Setter for the maximum allowed shift value determined by the smallest line
// /////////////////////////////////////////////////////////////////////////////
int Model :: set_dshift_max ()
{
    // Allocate memory
    dshift_max.resize(parameters->npoints());

    // For all points
    threaded_for(o, parameters->npoints(),
    {
        dshift_max[o] = std::numeric_limits<Real>::max();

        for (const LineProducingSpecies &lspec : lines.lineProducingSpecies)
        {
            const Real inverse_mass   = lspec.linedata.inverse_mass;
            const Real new_dshift_max = parameters->max_width_fraction
                                        * thermodynamics.profile_width (inverse_mass, o);

            if (dshift_max[o] > new_dshift_max)
            {
                dshift_max[o] = new_dshift_max;
            }
        }
    })

    return (0);
}
