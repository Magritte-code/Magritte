#include "parameters.hpp"
#include "tools/types.hpp"


#define READ_NUMBER(type, x)                                                               \
    type x##_copy;                                                                         \
    if (io.read_number ("."#x, x##_copy) == 0) {set_##x (x##_copy);}                       \
    else                                       {cout << "Failed to read "#x"!" << endl;}

#define READ_BOOL(type, x)                                                                 \
    type x##_copy;                                                                         \
    if (io.read_bool ("."#x, x##_copy) == 0) {set_##x (x##_copy);}                         \
    else                                     {cout << "Failed to read "#x"!" << endl;}

#define WRITE_NUMBER(type, x)                                                              \
    try         {io.write_number ("."#x, get_##x());}                                      \
    catch (...) {cout << "Failed write "#x"!" << endl;}

#define WRITE_BOOL(type, x)                                                                \
    try         {io.write_bool ("."#x, get_##x());}                                        \
    catch (...) {cout << "Failed write "#x"!" << endl;}


void Parameters :: read (const Io &io)
{
    cout << "Reading parameters..." << endl;

    READ_NUMBER (Size, npoints   );
    READ_NUMBER (Size, totnnbs   );
    READ_NUMBER (Size, nrays     );
    READ_NUMBER (Size, hnrays    );
    READ_NUMBER (Size, nrays_red );
    READ_NUMBER (Size, order_min );
    READ_NUMBER (Size, order_max );
    READ_NUMBER (Size, nboundary );
    READ_NUMBER (Size, nfreqs    );
    READ_NUMBER (Size, nspecs    );
    READ_NUMBER (Size, nlspecs   );
    READ_NUMBER (Size, nlines    );
    READ_NUMBER (Size, nquads    );

    READ_NUMBER (Real, pop_prec);

    READ_BOOL (bool, use_scattering      );
    READ_BOOL (bool, use_Ng_acceleration );
    READ_BOOL (bool, spherical_symmetry  );
    READ_BOOL (bool, adaptive_ray_tracing);

//
//
//    if   (io.read_number(".npoints",    n) == 0) {set_npoints   (n);}
//    else                                         {cout << "Failed read npoints!"    << endl;}
//    if   (io.read_number(".totnnbs",    n) == 0) {set_totnnbs   (n);}
//    else                                         {cout << "Failed read totnnbs!"    << endl;}
//    if   (io.read_number(".nrays",      n) == 0) {set_nrays     (n);}
//    else                                         {cout << "Failed read nrays!"      << endl;}
//    if   (io.read_number(".nrays_red",  n) == 0) {set_nrays_red (n);}
//    else                                         {cout << "Failed read nrays_red!"  << endl;}
//    if   (io.read_number(".order_min",  n) == 0) {set_order_min (n);}
//    else                                         {cout << "Failed read order_min!"  << endl;}
//    if   (io.read_number(".order_max",  n) == 0) {set_order_max (n);}
//    else                                         {cout << "Failed read order_max!"  << endl;}
//    if   (io.read_number(".nboundary",  n) == 0) {set_nboundary (n);}
//    else                                         {cout << "Failed read nboundary!"  << endl;}
//    if   (io.read_number(".nfreqs",     n) == 0) {set_nfreqs    (n);}
//    else                                         {cout << "Failed read nfreqs!"     << endl;}
//    if   (io.read_number(".nfreqs_red", n) == 0) {set_nfreqs_red(n);}
//    else                                         {cout << "Failed read nfreqs_red!" << endl;}
//    if   (io.read_number(".nspecs",     n) == 0) {set_nspecs    (n);}
//    else                                         {cout << "Failed read nspecs!"     << endl;}
//    if   (io.read_number(".nlspecs",    n) == 0) {set_nlspecs   (n);}
//    else                                         {cout << "Failed read nlspecs!"    << endl;}
//    if   (io.read_number(".nlines",     n) == 0) {set_nlines    (n);}
//    else                                         {cout << "Failed read nlines!"     << endl;}
//    if   (io.read_number(".nquads",     n) == 0) {set_nquads    (n);}
//    else                                         {cout << "Failed read nquads!"     << endl;}
//
//
//    double d;
//
//    if   (io.read_number(".pop_prec", d) == 0) {set_pop_prec(d);}
//    else                                       {cout << "Failed read pop_prec!" << endl;}
//
//
//    bool b;
//
//    if   (io.read_bool(".use_scattering",       b) == 0) {set_use_scattering       (b);}
//    else                                                 {cout << "Failed read use_scattering!"       << endl;}
//    if   (io.read_bool(".use_Ng_acceleration",  b) == 0) {set_use_Ng_acceleration  (b);}
//    else                                                 {cout << "Failed read use_Ng_acceleration!"  << endl;}
//    if   (io.read_bool(".spherical_symmetry",   b) == 0) {set_spherical_symmetry   (b);}
//    else                                                 {cout << "Failed read spherical_symmetry!"   << endl;}
//    if   (io.read_bool(".adaptive_ray_tracing", b) == 0) {set_adaptive_ray_tracing (b);}
//    else                                                 {cout << "Failed read adaptive_ray_tracing!" << endl;}
}




void Parameters :: write (const Io &io) const
{
    cout << "Writing parameters..." << endl;

//    try         {io.write_number (".npoints",    npoints    () );}
//    catch (...) {cout << "Failed write npoints!"         << endl;}
//    try         {io.write_number (".totnnbs",    totnnbs    () );}
//    catch (...) {cout << "Failed write totnnbs!"         << endl;}
//    try         {io.write_number (".nrays",      nrays      () );}
//    catch (...) {cout << "Failed write nrays!"           << endl;}
////    try         {io.write_number (".nrays_red",  nrays_red  () );}   // nrays_red depends on nprocs,
////    catch (...) {cout << "Failed write nrays_red!"       << endl;}   // so should be allowed to vary
//    try         {io.write_number (".order_min",  order_min  () );}
//    catch (...) {cout << "Failed write order_min!"       << endl;}
//    try         {io.write_number (".order_max",  order_max  () );}
//    catch (...) {cout << "Failed write order_max!"       << endl;}
//    try         {io.write_number (".nboundary",  nboundary  () );}
//    catch (...) {cout << "Failed write nboundary!"       << endl;}
//    try         {io.write_number (".nfreqs",     nfreqs    () );}
//    catch (...) {cout << "Failed write nfreqs!"          << endl;}
//    try         {io.write_number (".nfreqs_red", nfreqs_red () );}
//    catch (...) {cout << "Failed write nfreqs_red!"      << endl;}
//    try         {io.write_number (".nspecs",     nspecs     () );}
//    catch (...) {cout << "Failed write nspecs!"          << endl;}
//    try         {io.write_number (".nlspecs",    nlspecs    () );}
//    catch (...) {cout << "Failed write nlspecs!"         << endl;}
//    try         {io.write_number (".nlines",     nlines     () );}
//    catch (...) {cout << "Failed write nlines!"          << endl;}
//    try         {io.write_number (".nquads",     nquads     () );}
//    catch (...) {cout << "Failed write nquads!"          << endl;}

//    try         {io.write_number (".pop_prec", pop_prec() );}
//    catch (...) {cout << "Failed write pop_prec!"    << endl;}

//    try         {io.write_bool (".use_scattering",        use_scattering       () );}
//    catch (...) {cout << "Failed write use_scattering!"                     << endl;}
//    try         {io.write_bool (".use_Ng_acceleration",   use_Ng_acceleration  () );}
//    catch (...) {cout << "Failed write use_Ng_acceleration!"                << endl;}
//    try         {io.write_bool (".spherical_symmetry",    spherical_symmetry   () );}
//    catch (...) {cout << "Failed write spherical_symmetry!"                 << endl;}
//    try         {io.write_bool (".adaptive_ray_tracing",  adaptive_ray_tracing () );}
//    catch (...) {cout << "Failed write adaptive_ray_tracing!"               << endl;}

    WRITE_NUMBER (Size, npoints   );
    WRITE_NUMBER (Size, totnnbs   );
    WRITE_NUMBER (Size, nrays     );
    WRITE_NUMBER (Size, hnrays    );
    WRITE_NUMBER (Size, nrays_red );
    WRITE_NUMBER (Size, order_min );
    WRITE_NUMBER (Size, order_max );
    WRITE_NUMBER (Size, nboundary );
    WRITE_NUMBER (Size, nfreqs    );
    WRITE_NUMBER (Size, nspecs    );
    WRITE_NUMBER (Size, nlspecs   );
    WRITE_NUMBER (Size, nlines    );
    WRITE_NUMBER (Size, nquads    );

    WRITE_NUMBER (Real, pop_prec);

    WRITE_BOOL (bool, use_scattering      );
    WRITE_BOOL (bool, use_Ng_acceleration );
    WRITE_BOOL (bool, spherical_symmetry  );
    WRITE_BOOL (bool, adaptive_ray_tracing);
}
