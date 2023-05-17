#include "species.hpp"


const string prefix = "chemistry/species/";


///  read: read the input into the data structure
///    @param[in] io: io object
/////////////////////////////////////////////////
void Species :: read (const Io& io)
{
    cout << "Reading species..." << endl;

    parameters->set_nspecs  (io.get_length (prefix+"species"  ));
    parameters->set_npoints (io.get_length (prefix+"abundance"));

    abundance_init.resize (parameters->npoints());
    abundance     .resize (parameters->npoints());

    for (Size p = 0; p < parameters->npoints(); p++)
    {
      abundance_init[p].resize (parameters->nspecs());
      abundance     [p].resize (parameters->nspecs());
    }

    // Read the abundaces of each species in each cell
    io.read_array (prefix+"abundance", abundance);
    // Also read the species names
    io.read_list  (prefix+"species",   symbol   );

    // Set initial abundances
    abundance_init = abundance;
}


///  write: write out the data structure
///  @param[in] io: io object
////////////////////////////////////////
void Species :: write (const Io &io) const
{
    cout << "Writing species..." << endl;

    Long1 dummy (parameters->nspecs(), 0);

    io.write_list  (prefix+"species",   symbol   );
    io.write_array (prefix+"abundance", abundance);
}
