#include "frequencies.hpp"
#include "Tools/constants.hpp"
#include "Tools/types.hpp"


const string prefix = "radiation/frequencies/";


///  Reader for the Frequencies data
///    @param[in] io : io object to read with
/////////////////////////////////////////////
void Frequencies :: read (const Io& io)
{
    cout << "Reading frequencies..." << endl;

    cout << "ncells = " << ncells << endl;
    cout << "nlines = " << nlines << endl;
    cout << "nquads = " << nquads << endl;

    // Count line frequencies
    nfreqs = nlines * nquads;

    // Add extra frequency bins around lines to get nicer spectrum
    //nfreqs += nlines * 2 * nbins;

    // Add ncont bins background
    //nfreqs += ncont;

    // Ensure that nfreq is a multiple of n_simd_lanes
    nfreqs_red = reduced (nfreqs);
    nfreqs     = nfreqs_red * n_simd_lanes;

    nu.resize (npoints);

    for (Size p = 0; p < npoints; p++)
    {
        nu[p].resize (nfreqs_red);
    }

    appears_in_line_integral.resize (nfreqs);
    corresponding_l_for_spec.resize (nfreqs);
    corresponding_k_for_tran.resize (nfreqs);
    corresponding_z_for_line.resize (nfreqs);

    // frequencies.nu has to be initialized (for unused entries)
    threaded_for (p, npoints)
    {
        for (Size f = 0; f < nfreqs_red; f++)
        {
            nu[p][f] = 0.0;
        }
    }
}


///  Writer for the Frequencies data
///    @param[in] io : io object to write with
/////////////////////.////////////////////////
void Frequencies :: write (const Io &io) const
{
    cout << "Writing frequencies..." << endl;

    io.write_array (prefix+"nu", nu);
}
