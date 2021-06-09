#include "collisionPartner.hpp"
#include "tools/types.hpp"


const string prefix = "lines/lineProducingSpecies_";


///  read: read in collision partner data
///    @param[in] io: io object
/////////////////////////////////////////
void CollisionPartner :: read (const Io& io, const Size l, const Size c)
{
    cout << "Reading collisionPartner..." << endl;

    const string prefix_lc = prefix + std::to_string (l) + "/linedata"
                             + "/collisionPartner_" + std::to_string (c) + "/";

    io.read_number (prefix_lc+".num_col_partner", num_col_partner);
    io.read_word   (prefix_lc+".orth_or_para_H2", orth_or_para_H2);

    io.read_length (prefix_lc+"tmp",  ntmp);
    io.read_length (prefix_lc+"icol", ncol);

    tmp.resize (ntmp);
    io.read_list (prefix_lc+"tmp", tmp);

    icol.resize (ncol);
    jcol.resize (ncol);
    io.read_list (prefix_lc+"icol", icol);
    io.read_list (prefix_lc+"jcol", jcol);

    Ce.resize (ntmp);
    Cd.resize (ntmp);

    for (Size t = 0; t < ntmp; t++)
    {
        Ce[t].resize (ncol);
        Cd[t].resize (ncol);
    }

    io.read_array (prefix_lc+"Ce", Ce);
    io.read_array (prefix_lc+"Cd", Cd);

    Ce_intpld.resize (ncol);
    Cd_intpld.resize (ncol);
}


///  write: read in collision partner data
///    @param[in] io: io object
//////////////////////////////////////////
void CollisionPartner :: write (const Io& io, const Size l, const Size c) const
{
    cout << "Writing collisionPartner..." << endl;
    cout << "(l, c) = " << l << ", " << c << endl;

    const string prefix_lc = prefix + std::to_string (l) + "/linedata"
                             + "/collisionPartner_" + std::to_string (c) + "/";

    io.write_number (prefix_lc+".num_col_partner", num_col_partner);
    io.write_word   (prefix_lc+".orth_or_para_H2", orth_or_para_H2);

    io.write_list (prefix_lc+"icol", icol);
    io.write_list (prefix_lc+"jcol", jcol);

    io.write_list (prefix_lc+"tmp", tmp);

    io.write_array (prefix_lc+"Ce", Ce);
    io.write_array (prefix_lc+"Cd", Cd);
}
