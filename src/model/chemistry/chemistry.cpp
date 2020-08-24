#include "chemistry.hpp"


///  Reader for model data
///    @param[in] io: io data object
////////////////////////////////////
void Chemistry :: read (const Io& io)
{
  cout << "Reading chemistry..." << endl;

  species.read (io);
}


///  Writer for model data
///    @param[in] io: io data object
////////////////////////////////////
void Chemistry :: write (const Io& io) const
{
  cout << "Writing chemistry..." << endl;

  species.write (io);
}
