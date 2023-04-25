#include "parameters.hpp"
#include "tools/types.hpp"


#define READ_SETONCE_NUMBER(type, x, essential)                                            \
    type x##_copy;                                                               \
    if (io.read_number ("."#x, x##_copy) == 0) {set_##x (x##_copy);}                       \
    else                                       {if(essential){cout << "Failed to read "#x"!" << endl;}}

#define READ_SETONCE_BOOL(type, x, essential)                                              \
    type x##_copy;                                                                         \
    if (io.read_bool ("."#x, x##_copy) == 0) {set_##x (x##_copy);}                         \
    else                                     {if(essential){cout << "Failed to read "#x"!" << endl;}}

#define READ_SETONCE_WORD(type, x, essential)                                              \
    type x##_copy;                                                                         \
    if (io.read_word ("."#x, x##_copy) == 0) {set_##x (x##_copy);}                         \
    else                                     {if(essential){cout << "Failed to read "#x"!" << endl;}}

#define READ_NUMBER(type, x)                                                               \
    type x##_copy;                                                                         \
    if (io.read_number ("."#x, x##_copy) == 0) {x = x##_copy;}                             \
    else                                       {cout << "Failed to read "#x"!" << endl;}

#define READ_BOOL(type, x)                                                                 \
    type x##_copy;                                                                         \
    if (io.read_bool ("."#x, x##_copy) == 0) {x = x##_copy;}                               \
    else                                     {cout << "Failed to read "#x"!" << endl;}

#define WRITE_SETONCE_NUMBER(type, x, essential)                                           \
    try         {io.write_number ("."#x, x());}                                            \
    catch (...) {if(essential){cout << "Failed to write parameter "#x"!" << endl;}}

#define WRITE_SETONCE_BOOL(type, x, essential)                                             \
    try         {io.write_bool ("."#x, x());}                                              \
    catch (...) {if(essential){cout << "Failed to write parameter "#x"!" << endl;}}

#define WRITE_SETONCE_WORD(type, x, essential)                                             \
    try         {io.write_word ("."#x, x());}                                              \
    catch (...) {if(essential){cout << "Failed to write parameter "#x"!" << endl;}}

#define WRITE_NUMBER(type, x)                                                              \
    try         {io.write_number ("."#x, x);}                                              \
    catch (...) {cout << "Failed to write parameter "#x"!" << endl;}

#define WRITE_BOOL(type, x)                                                                \
    try         {io.write_bool ("."#x, x);}                                                \
    catch (...) {cout << "Failed to write parameter "#x"!" << endl;}


void Parameters :: read (const Io &io)
{
    cout << "Reading parameters..." << endl;

    READ_SETONCE_NUMBER (Size, dimension, true);
    READ_SETONCE_NUMBER (Size, npoints, true);
    READ_SETONCE_NUMBER (Size, nrays, true);
    READ_SETONCE_NUMBER (Size, hnrays, false);
    READ_SETONCE_NUMBER (Size, nboundary, true);
    READ_SETONCE_NUMBER (Size, nfreqs, false);
    READ_SETONCE_NUMBER (Size, nspecs, true);
    READ_SETONCE_NUMBER (Size, nlspecs, true);
    READ_SETONCE_NUMBER (Size, nlines, false);
    READ_SETONCE_NUMBER (Size, nquads, true);

    READ_SETONCE_BOOL (bool, use_scattering, false);
    READ_SETONCE_BOOL (bool, spherical_symmetry, true);
    READ_SETONCE_BOOL (bool, adaptive_ray_tracing, false);

    READ_NUMBER (Size, n_off_diag);

    READ_NUMBER (Real, max_width_fraction         );
    READ_NUMBER (Real, convergence_fraction       );
    READ_NUMBER (Real, min_rel_pop_for_convergence);
    READ_NUMBER (Real, pop_prec                   );
    READ_NUMBER (Real, min_opacity                );
    READ_NUMBER (Real, Ng_acceleration_mem_limit);

    READ_BOOL (bool, store_intensities     );
    READ_BOOL (bool, one_line_approximation);
    READ_BOOL (bool, use_adaptive_Ng_acceleration);
}




void Parameters :: write (const Io &io) const
{
    cout << "Writing parameters..." << endl;

    WRITE_SETONCE_WORD (string, version, true);

    WRITE_SETONCE_NUMBER (Size, dimension, true);
    WRITE_SETONCE_NUMBER (Size, npoints, true);
    WRITE_SETONCE_NUMBER (Size, nrays, true);
    WRITE_SETONCE_NUMBER (Size, hnrays, false);
    WRITE_SETONCE_NUMBER (Size, nboundary, true);
    WRITE_SETONCE_NUMBER (Size, nfreqs, false);
    WRITE_SETONCE_NUMBER (Size, nspecs, true);
    WRITE_SETONCE_NUMBER (Size, nlspecs, true);
    WRITE_SETONCE_NUMBER (Size, nlines, false);
    WRITE_SETONCE_NUMBER (Size, nquads, true);


    WRITE_SETONCE_BOOL (bool, use_scattering, false);
    WRITE_SETONCE_BOOL (bool, spherical_symmetry, true);
    WRITE_SETONCE_BOOL (bool, adaptive_ray_tracing, false);

    WRITE_NUMBER (Size, n_off_diag);

    WRITE_NUMBER (Real, max_width_fraction         );
    WRITE_NUMBER (Real, convergence_fraction       );
    WRITE_NUMBER (Real, min_rel_pop_for_convergence);
    WRITE_NUMBER (Real, pop_prec                   );
    WRITE_NUMBER (Real, min_opacity                );
    WRITE_NUMBER (Real, Ng_acceleration_mem_limit);

    WRITE_BOOL (bool, store_intensities     );
    WRITE_BOOL (bool, one_line_approximation);
    WRITE_BOOL (bool, use_adaptive_Ng_acceleration);
}
