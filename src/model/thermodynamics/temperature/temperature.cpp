#include "temperature.hpp"

const string prefix = "thermodynamics/temperature/";

///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Temperature ::read(const Io& io) {
    cout << "Reading temperature..." << endl;

    parameters->set_npoints(io.get_length(prefix + "gas"));

    gas.resize(parameters->npoints());

    io.read_list(prefix + "gas", gas);
}

///  write: write out data structure
///    @param[in] io: io object
/////////////////////////////////
void Temperature ::write(const Io& io) const {
    cout << "Writing temperature..." << endl;

    io.write_list(prefix + "gas", gas);
}
