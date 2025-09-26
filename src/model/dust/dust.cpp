#include "dust.hpp"

/// Reader for dust data
///    @param[in] io: io data object
void Dust::read(const Io& io) {
    cout << "Reading dust..." << endl;

    // by default initialized to zero
    n_dust_frequencies = 0;
    // Note: dust is optional to include in the model, so first check if data is present
    int err = io.read_number("dust/.n_dust_frequencies", n_dust_frequencies);
    if (err == 0 && n_dust_frequencies > 0) {
        dust_temperature.resize(parameters->npoints());
        dust_frequencies.resize(n_dust_frequencies);
        dust_opacities.resize(parameters->npoints(), n_dust_frequencies);

        io.read_list("dust/dust_temperature", dust_temperature);
        io.read_list("dust/dust_frequencies", dust_frequencies);
        io.read_list("dust/dust_opacities", dust_opacities);

        cout << "  Read dust/continuum opacities at " << n_dust_frequencies << " frequencies."
             << endl;
    } else {
        // in this case, n_dust_frequencies remains zero, so when calculating dust properties, this
        // must be checked to skip that step in case of no dust present in the model.
        cout << "No dust/continuum data found." << endl;
    }
}

/// Writer for dust data
///    @param[in] io: io data object
void Dust::write(const Io& io) const {
    cout << "Writing dust..." << endl;

    // Note: data will be written either way, even if it all zeros
    io.write_list("dust/dust_temperature", dust_temperature);
    io.write_list("dust/dust_frequencies", dust_frequencies);
    io.write_list("dust/dust_opacities", dust_opacities);

    Size temp_n_dust_frequencies = dust_frequencies.size();

    io.write_number("dust/.n_dust_frequencies", temp_n_dust_frequencies);
}
