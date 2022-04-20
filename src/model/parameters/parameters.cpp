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

    READ_NUMBER (Size, dimension );
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
    READ_BOOL (bool, store_intensities   );
    READ_BOOL (bool, use_Ng_acceleration );
    READ_BOOL (bool, spherical_symmetry  );
    READ_BOOL (bool, adaptive_ray_tracing);
}




void Parameters :: write (const Io &io) const
{
    cout << "Writing parameters..." << endl;

    WRITE_NUMBER (Size, dimension );
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
    WRITE_BOOL (bool, store_intensities   );
    WRITE_BOOL (bool, use_Ng_acceleration );
    WRITE_BOOL (bool, spherical_symmetry  );
    WRITE_BOOL (bool, adaptive_ray_tracing);
}
