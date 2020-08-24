#include "quadrature.hpp"


const string prefix = "Lines/LineProducingSpecies_";


///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Quadrature :: read (const Io &io, const Size l)
{
    cout << "Reading quadrature..." << endl;

    const string prefix_l = prefix + std::to_string (l) + "/quadrature/";

    io.read_length (prefix_l+"weights", nquads);

    weights.resize (nquads);
    roots  .resize (nquads);

    io.read_list (prefix_l+"weights", weights);
    io.read_list (prefix_l+"roots",   roots  );
}


///  write: write out data structure
///    @param[in] io: io object
////////////////////////////////////
void Quadrature :: write (const Io &io, const Size l) const
{
    cout << "Writing quadrature..." << endl;

    const string prefix_l = prefix + std::to_string (l) + "/quadrature/";

    io.write_list (prefix_l+"weights", weights);
    io.write_list (prefix_l+"roots",   roots  );
}
