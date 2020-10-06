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

    cout << "npoints = " << parameters.npoints() << endl;
    cout << "nlines  = " << parameters.nlines()  << endl;
    cout << "nquads  = " << parameters.nquads()  << endl;

    // Count line frequencies
    parameters.set_nfreqs (parameters.nlines() * parameters.nquads());

    // Add extra frequency bins around lines to get nicer spectrum
    //nfreqs += nlines * 2 * nbins;

    // Add ncont bins background
    //nfreqs += ncont;

    nu.resize (parameters.npoints(), parameters.nfreqs());

    appears_in_line_integral.resize (parameters.nfreqs());
    corresponding_l_for_spec.resize (parameters.nfreqs());
    corresponding_k_for_tran.resize (parameters.nfreqs());
    corresponding_z_for_line.resize (parameters.nfreqs());

    // frequencies.nu has to be initialized (for unused entries)
    threaded_for (p, parameters.npoints(),
    {
        for (Size f = 0; f < parameters.nfreqs(); f++)
        {
            nu(p,f) = 0.0;
        }
    })
}


///  Writer for the Frequencies data
///    @param[in] io : io object to write with
/////////////////////.////////////////////////
void Frequencies :: write (const Io &io) const
{
    cout << "Writing frequencies..." << endl;

    io.write_array (prefix+"nu", nu);
}
