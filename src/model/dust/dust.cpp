#include "dust.hpp"
#include "tools/constants.hpp"
#include "tools/types.hpp"


const string prefix = "dust/";


///  read: read in dust data
///    @param[in] io: io object
///////////////////////////////
void Dust :: read (const Io& io)
{
    cout << "Reading dust data..." << endl;

    Size num_f = 0;
    io.read_number (prefix+".num_f", num_f);

    cout << "read num_f " << num_f << endl;

    freqs.resize (num_f);
    kappa.resize (num_f);

    io.read_list (prefix+"freqs", freqs);
    io.read_list (prefix+"kappa", kappa);

    density.resize (parameters.npoints());

    io.read_list (prefix+"density", density);
}


///  write: write out data structure
///    @param[in] io: io object
////////////////////////////////////
void Dust :: write (const Io& io) const
{
    cout << "Writing dust data..." << endl;

    const Size num_f = freqs.size();
    io.write_number (prefix+".num_f", num_f);

    io.write_list (prefix+"freqs", freqs);
    io.write_list (prefix+"kappa", kappa);

    io.write_list (prefix+"density", density);
}