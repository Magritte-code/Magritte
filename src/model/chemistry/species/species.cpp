#include "species.hpp"


const string prefix = "chemistry/species/";


///  read: read the input into the data structure
///    @param[in] io: io object
/////////////////////////////////////////////////
void Species :: read (const Io& io)
{
    cout << "Reading species..." << endl;

    parameters.set_nspecs  (io.get_length (prefix+"species"  ));
    cout << "nspecs  = " <<  parameters.nspecs () << endl;

    parameters.set_npoints (io.get_length (prefix+"abundance"));
    cout << "npoints = " <<  parameters.npoints() << endl;

    abundance_init.resize (parameters.npoints());
    cout << "Here? 1" << endl;
    abundance     .resize (parameters.npoints());

    cout << "Here? 2" << endl;

    for (Size p = 0; p < parameters.npoints(); p++)
    {
      abundance_init[p].resize (parameters.nspecs());
      abundance     [p].resize (parameters.nspecs());
    }

    cout << "Here? 3" << endl;

    // Read the abundaces of each species in each cell
    io.read_array (prefix+"abundance", abundance);

    cout << "Here? 4" << endl;

    // Set initial abundances
    abundance_init = abundance;

    cout << "Here? 5" << endl;
}


///  write: write out the data structure
///  @param[in] io: io object
////////////////////////////////////////
void Species :: write (const Io &io) const
{
    cout << "Writing species..." << endl;

    Long1 dummy (parameters.nspecs(), 0);

    io.write_list  (prefix+"species",   dummy    );
    io.write_array (prefix+"abundance", abundance);
}
