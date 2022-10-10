#include "cooling.hpp"
const string prefix = "cooling/";

// should contain the setup for the cooling storage
// for all points (and all freqs?), some cooling rate should be computed
// TODO: figure out which to use
// 1) store for all freqs: no reasonable formulation for the collisional approach
// 2) store sum/integral: not able to pinpoint what contributes most to cooling...
// currently, option (2) seems best.


///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Cooling :: read (const Io& io)
{
    cooling_rate.resize(parameters->npoints());
    io.read_list(prefix+"cooling_rate", cooling_rate);

}

// // Initialises the structure for computing the cooling rates
// void Cooling :: setup (Model& model)
// {
//     cooling_rate.resize(model.parameters->npoints());
// }


///  Writer for the Cooling data
///    @param[in] io : io object to write with
/////////////////////.////////////////////////
void Cooling :: write (const Io &io) const
{
    cout << "Writing cooling rates..." << endl;

    io.write_list(prefix+"cooling_rate", cooling_rate);
}
