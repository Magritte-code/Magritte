#include "turbulence.hpp"

const string prefix = "thermodynamics/turbulence/";

///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Turbulence ::read(const Io& io) {
    cout << "Reading turbulence..." << endl;

    parameters->set_npoints(io.get_length(prefix + "vturb2"));

    vturb2.resize(parameters->npoints());

    io.read_list(prefix + "vturb2", vturb2);
}

///  write: write out data structure
///    @param[in] io: io object
/////////////////////////////////
void Turbulence ::write(const Io& io) const {
    cout << "Writing turbulence..." << endl;

    io.write_list(prefix + "vturb2", vturb2);
}
