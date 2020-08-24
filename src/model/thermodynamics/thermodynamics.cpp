#include "thermodynamics.hpp"


///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Thermodynamics :: read (const Io& io)
{
    cout << "Reading thermodynamics..." << endl;

    temperature.read (io);
    turbulence .read (io);
}


///  write: write out data structure
///    @param[in] io: io object
////////////////////////////////////
void Thermodynamics :: write (const Io& io) const
{
    cout << "Writing thermodynamics..." << endl;

    temperature.write (io);
    turbulence .write (io);
}
