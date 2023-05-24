#include "linedata.hpp"

#include "tools/constants.hpp"
#include "tools/types.hpp"

const string prefix = "lines/lineProducingSpecies_";

///  read: read in line data
///    @param[in] io: io object
///    @param[in] l: nr of line producing species
/////////////////////////////////////////////////
void Linedata ::read(const Io& io, const Size l) {
    cout << "Reading linedata..." << endl;

    const string prefix_l = prefix + std::to_string(l) + "/linedata/";

    io.read_number(prefix_l + ".num", num);

    cout << "read num " << num << endl;

    io.read_word(prefix_l + ".sym", sym);

    cout << "read sym " << sym << endl;

    io.read_number(prefix_l + ".inverse_mass", inverse_mass);

    io.read_number(prefix_l + ".nlev", nlev);
    io.read_number(prefix_l + ".nrad", nrad);

    cout << "nlev = " << nlev << endl;
    cout << "nrad = " << nrad << endl;

    irad.resize(nrad);
    jrad.resize(nrad);

    io.read_list(prefix_l + "irad", irad);
    io.read_list(prefix_l + "jrad", jrad);

    energy.resize(nlev);
    weight.resize(nlev);

    io.read_list(prefix_l + "energy", energy);
    io.read_list(prefix_l + "weight", weight);

    frequency.resize(nrad);

    io.read_list(prefix_l + "frequency", frequency);

    A.resize(nrad);
    Bs.resize(nrad);
    Ba.resize(nrad);

    io.read_list(prefix_l + "A", A);
    io.read_list(prefix_l + "Bs", Bs);
    io.read_list(prefix_l + "Ba", Ba);

    // Get ncolpar
    io.read_length(prefix_l + "collisionPartner_", ncolpar);

    colpar.resize(ncolpar);

    for (Size c = 0; c < ncolpar; c++) {
        colpar[c].read(io, l, c);
    }

    ncol_tot = 0;

    for (Size c = 0; c < ncolpar; c++) {
        ncol_tot += colpar[c].ncol;
    }
}

///  write: write out data structure
///    @param[in] io: io object
///    @param[in] l: nr of line producing species
/////////////////////////////////////////////////
void Linedata ::write(const Io& io, const Size l) const {
    cout << "Writing linedata..." << endl;

    const string prefix_l = prefix + std::to_string(l) + "/linedata/";

    io.write_number(prefix_l + ".num", num);
    io.write_word(prefix_l + ".sym", sym);

    io.write_number(prefix_l + ".inverse_mass", inverse_mass);

    io.write_number(prefix_l + ".nlev", nlev);
    io.write_number(prefix_l + ".nrad", nrad);

    io.write_list(prefix_l + "irad", irad);
    io.write_list(prefix_l + "jrad", jrad);

    io.write_list(prefix_l + "energy", energy);
    io.write_list(prefix_l + "weight", weight);

    io.write_list(prefix_l + "frequency", frequency);

    io.write_list(prefix_l + "A", A);
    io.write_list(prefix_l + "Bs", Bs);
    io.write_list(prefix_l + "Ba", Ba);

    cout << "ncolpoar = " << ncolpar << endl;

    for (Size c = 0; c < ncolpar; c++) {
        cout << "--- colpoar = " << c << endl;
        colpar[c].write(io, l, c);
    }
}
