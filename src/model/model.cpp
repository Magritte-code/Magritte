#include "model.hpp"
#include "tools/heapsort.hpp"
#include "paracabs.hpp"


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
    cout << "  nfreqs_red = " << parameters.nfreqs_red () << endl;
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
