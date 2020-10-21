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
int Model :: compute_radiation_field_0th_short_characteristics ()
{
    cout << "Computing radiation field..." << endl;

    const Size length_max = 4*parameters.npoints() + 1;
    const Size  width_max =   parameters.nfreqs ();

    Solver solver (length_max, width_max, parameters.n_off_diag);
    solver.solve_0th_order_short_charateristics  (*this);

    return (0);
}


///  Computer for the radiation field
/////////////////////////////////////
int Model :: compute_radiation_field_2nd_order_Feautrier ()
{
    cout << "Computing radiation field..." << endl;

    const Size length_max = 4*parameters.npoints() + 1;
    const Size  width_max =   parameters.nfreqs ();


    cout << "npoints = " << parameters.npoints() << endl;
    cout << "nfreqs  = " << parameters.nfreqs () << endl;

    cout << "l_max = " << length_max << endl;
    cout << "w_max = " <<  width_max << endl;

    Solver solver (length_max, width_max, parameters.n_off_diag);
    solver.solve_2nd_order_Feautrier (*this);

    return (0);
}


///  Compute the effective mean intensity in a line
///////////////////////////////////////////////////
int Model :: compute_Jeff ()
{
    for (LineProducingSpecies &lspec : lines.lineProducingSpecies)
    {
//        Lambda = MatrixXd::Zero (lspec.population.size(), lspec.population.size());

        threaded_for (p, parameters.npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                const Long1 freq_nrs = lspec.nr_line[p][k];

                // Initialize values
                lspec.Jlin[p][k] = 0.0;

                // Integrate over the line
                for (Size z = 0; z < parameters.nquads(); z++)
                {
                    lspec.Jlin[p][k] += lspec.quadrature.weights[z] * radiation.J(p, freq_nrs[z]);
                }


                double diff = 0.0;

                // Collect the approximated part
//                for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
//                {
//                    const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);
//
//                    diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
//                }

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
    lines.iteration_using_statistical_equilibrium (
            chemistry.species.abundance,
            thermodynamics.temperature.gas,
            parameters.pop_prec()                 );

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
        // const Io   &io,
        const bool  use_Ng_acceleration,
        const long  max_niterations     )
{
    // Check spectral discretisation setting
    if (spectralDiscretisation != SD_Lines)
    {
        throw std::runtime_error ("Spectral discretisation was not set for Lines!");
    }

    // Write out initial level populations
    //for (int l = 0; l < parameters.nlspecs(); l++)
    //{
    //    const string tag = "_rank_" + str_MPI_comm_rank() + "_iteration_0";

    //    lines.lineProducingSpecies[l].write_populations (io, l, tag);
    //}

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

        //if (use_Ng_acceleration && (iteration_normal == 4))
        //{
        //    lines.iteration_using_Ng_acceleration (
        //        parameters.pop_prec()             );
        //
        //    iteration_normal = 0;
        //}
        // else
        {
            // logger.write ("Computing the radiation field...");
            cout << "Computing the radiation field..." << endl;

            compute_radiation_field_2nd_order_Feautrier ();
            compute_Jeff            ();

            lines.iteration_using_statistical_equilibrium (
                chemistry.species.abundance,
                thermodynamics.temperature.gas,
                parameters.pop_prec()                     );

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

            // const string tag = "_rank_" + str_MPI_comm_rank() + "_iteration_" + to_string (iteration);

            // lines.lineProducingSpecies[l].write_populations (io, l, tag);
        }
    } // end of while loop of iterations

    // Print convergence stats
    cout << "Converged after " << iteration << " iterations" << endl;

    return iteration;
}
