#include "species.hpp"


const string prefix = "chemistry/species/";


///  read: read the input into the data structure
///    @param[in] io: io object
/////////////////////////////////////////////////
void Species :: read (const Io& io)
{
    cout << "Reading species..." << endl;

    io.read_length (prefix+"species",   nspecs);
    io.read_length (prefix+"abundance", npoints);

    abundance_init.resize (npoints);
    abundance     .resize (npoints);

    for (Size p = 0; p < ncells; p++)
    {
      abundance_init[p].resize (nspecs);
      abundance     [p].resize (nspecs);
    }

    // Read the abundaces of each species in each cell
    io.read_array (prefix+"abundance", abundance);

    // Set initial abundances
    abundance_init = abundance;
}


///  write: write out the data structure
///  @param[in] io: io object
////////////////////////////////////////
void Species :: write (const Io &io) const
{
    cout << "Writing species..." << endl;

    Long1 dummy (nspecs, 0);

    io.write_list  (prefix+"species",   dummy    );
    io.write_array (prefix+"abundance", abundance);
}
