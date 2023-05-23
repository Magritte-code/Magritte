#include "frequencies.hpp"

#include "tools/constants.hpp"
#include "tools/types.hpp"

const string prefix = "radiation/frequencies/";

///  Reader for the Frequencies data
///    @param[in] io : io object to read with
/////////////////////////////////////////////
void Frequencies ::read(const Io& io) {
    cout << "Reading frequencies..." << endl;

    // Count line frequencies
    parameters->set_nfreqs(parameters->nlines() * parameters->nquads());

    resize_data(parameters->nfreqs());
}

/// Helper function for resize all data associated to
void Frequencies ::resize_data(const Size nfreqs) {
    nu.resize(parameters->npoints(), nfreqs);

    corresponding_line.resize(nfreqs);

    appears_in_line_integral.resize(nfreqs);
    corresponding_l_for_spec.resize(nfreqs);
    corresponding_k_for_tran.resize(nfreqs);
    corresponding_z_for_line.resize(nfreqs);

    // frequencies.nu has to be initialized (for unused entries)
    threaded_for (p, parameters->npoints(), {
        for (Size f = 0; f < nfreqs; f++) {
            nu(p, f) = 0.0;
        }
    })

        parameters->set_nfreqs(nfreqs); // is global variable, so may be reset to another value
}

///  Writer for the Frequencies data
///    @param[in] io : io object to write with
/////////////////////.////////////////////////
void Frequencies ::write(const Io& io) const {
    cout << "Writing frequencies..." << endl;

    io.write_array(prefix + "nu", nu);
}
