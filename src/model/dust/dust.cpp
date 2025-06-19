#include "dust.hpp"

/// Reader for dust data
///    @param[in] io: io data object
void Dust::read(const Io& io) {
    cout << "Reading dust..." << endl;

    // by default initialized to zero
    dust_opacities.resize(parameters->npoints(), parameters->nlines());
    dust_emissivities.resize(parameters->npoints(), parameters->nlines());

    // Note: dust is optional to include in the model, so first check if data is present
    Matrix<Real> temp_dust_opacities(parameters->npoints(), parameters->nlines());
    int err = io.read_array("dust/opacities", temp_dust_opacities);

    if (err == 0) {
        // TODO: check if the data can be read this way
        io.read_array("dust/opacities", dust_opacities);
        io.read_array("dust/emissivities", dust_emissivities);
    }
}

/// Writer for dust data
///    @param[in] io: io data object
void Dust::write(const Io& io) const {
    cout << "Writing dust..." << endl;

    // Note: data will be written either way, even if it all zeros
    io.write_array("dust/opacities", dust_opacities);
    io.write_array("dust/emissivities", dust_emissivities);
}