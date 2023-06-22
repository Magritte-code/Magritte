#include "quadrature.hpp"

const string prefix = "lines/lineProducingSpecies_";

///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Quadrature ::read(const Io& io, const Size l) {
    cout << "Reading quadrature..." << endl;

    const string prefix_l = prefix + std::to_string(l) + "/quadrature/";

    parameters->set_nquads(io.get_length(prefix_l + "weights"));

    weights.resize(parameters->nquads());
    roots.resize(parameters->nquads());

    io.read_list(prefix_l + "weights", weights);
    io.read_list(prefix_l + "roots", roots);
}

///  write: write out data structure
///    @param[in] io: io object
////////////////////////////////////
void Quadrature ::write(const Io& io, const Size l) const {
    cout << "Writing quadrature..." << endl;

    const string prefix_l = prefix + std::to_string(l) + "/quadrature/";

    io.write_list(prefix_l + "weights", weights);
    io.write_list(prefix_l + "roots", roots);
}
