#include "frequencies.hpp"
#include "tools/constants.hpp"
#include "tools/types.hpp"


const string prefix = "radiation/frequencies/";


///  Reader for the Frequencies data
///    @param[in] io : io object to read with
/////////////////////////////////////////////
void Frequencies :: read (const Io& io)
{
    cout << "Reading frequencies..." << endl;

    // Count line frequencies
    parameters->set_nfreqs (parameters->nlines() * parameters->nquads());

    resize_data(parameters->nfreqs());
}

/// Helper function for resize all data associated to
void Frequencies :: resize_data (const Size Nfreqs)
{
    nu.resize (parameters->npoints(), Nfreqs);

    corresponding_line.resize (Nfreqs);

    appears_in_line_integral.resize (Nfreqs);
    corresponding_l_for_spec.resize (Nfreqs);
    corresponding_k_for_tran.resize (Nfreqs);
    corresponding_z_for_line.resize (Nfreqs);

    // frequencies.nu has to be initialized (for unused entries)
    threaded_for (p, parameters->npoints(),
    {
        for (Size f = 0; f < Nfreqs; f++)
        {
            nu(p,f) = 0.0;
        }
    })

    N_IMAGE_FREQS = Nfreqs;
}


///  Writer for the Frequencies data
///    @param[in] io : io object to write with
/////////////////////.////////////////////////
void Frequencies :: write (const Io &io) const
{
    cout << "Writing frequencies..." << endl;

    io.write_array (prefix+"nu", nu);
}
