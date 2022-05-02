#pragma once


#include <limits>

#include "io/io.hpp"
#include "tools/setOnce.hpp"


///  Create for each parameter "x":
///    - a (private) local variable "x__",
///    - a (private) Meyers' singleton for isSet "x_isSet()",
///    - a (public) Meyers' singleton for value "x()",
///    - a (public) setter function "set_x".
///  Remember to ensure in Parameters copy constructor to
///  ensure that the local variables are up to date, i.e.
///  ensure that their values equal the global values.

#define CREATE_PARAMETER(type, x)                                                           \
    private:                                                                                \
        SetOnce<type> x##__;                     /* Local (SetOnce) variable           */   \
        accel inline bool& x##_isSet () const    /* Meyers' singleton for global isSet */   \
        {                                                                                   \
            static bool isSet = false;           /* Global isSet                       */   \
            return isSet;                        /* Return global isSet                */   \
        }                                                                                   \
    public:                                                                                 \
        accel inline type& x () const            /* Meyers' singleton for global value */   \
        {                                                                                   \
            static type x##_value;               /* Global value                       */   \
            if (!x##_isSet())                    /* If global value is not set         */   \
            {                                                                               \
                x##_value   = x##__.get();       /* Set global with local value        */   \
                x##_isSet() = true;              /* Mark global value as set           */   \
            }                                                                               \
            return x##_value;                    /* Return global value                */   \
        }                                                                                   \
        inline void set_##x (const type value)   /* Setter function                    */   \
        {                                                                                   \
            x##__.set (value);                   /* Set local value                    */   \
            x() = value;                         /* Set global value                   */   \
        }                                                                                   \
        inline type get_##x () const             /* Getter function                    */   \
        {                                                                                   \
            return x();                          /* Return copy of global value        */   \
        }

#define CONSTRUCT_PARAMETER(type, x)                                                        \
    x##__ = SetOnce <type> (x());

#define COPY_PARAMETER(x)                                                                   \
    if (x##_isSet()) x##__.set(x());



///  Parameters: secure structure for the model parameters
//////////////////////////////////////////////////////////
struct Parameters
{
    long   n_off_diag           = 0;
    double max_width_fraction   = 0.5;
    double convergence_fraction = 0.995;

    Size max_matrix_size=100*1024*1024;///< Max temp matrix size when solving for the level populations //By default hundred megabyte
    //Max extra memory usage during stat eq solving is then approx times 2 due to use of eigen triplets (1) + setting matrix from this (1).

    void read (const Io &io);
    void write(const Io &io) const;

    CREATE_PARAMETER (string, model_name);

    CREATE_PARAMETER (Size, dimension );
    CREATE_PARAMETER (Size, npoints   );
    CREATE_PARAMETER (Size, totnnbs   );
    CREATE_PARAMETER (Size, nrays     );
    CREATE_PARAMETER (Size, hnrays    );
    CREATE_PARAMETER (Size, order_min );
    CREATE_PARAMETER (Size, order_max );
    CREATE_PARAMETER (Size, nboundary );
    CREATE_PARAMETER (Size, nfreqs    );
    CREATE_PARAMETER (Size, nspecs    );
    CREATE_PARAMETER (Size, nlspecs   );
    CREATE_PARAMETER (Size, nlines    );
    CREATE_PARAMETER (Size, nquads    );

    CREATE_PARAMETER (Real, pop_prec);

    CREATE_PARAMETER (bool, use_scattering      );
    CREATE_PARAMETER (bool, store_intensities   );
    CREATE_PARAMETER (bool, use_Ng_acceleration );
    CREATE_PARAMETER (bool, spherical_symmetry  );
    CREATE_PARAMETER (bool, adaptive_ray_tracing);

    Parameters ()
    {
        CONSTRUCT_PARAMETER (string, model_name);

        CONSTRUCT_PARAMETER (Size, dimension );
        CONSTRUCT_PARAMETER (Size, npoints   );
        CONSTRUCT_PARAMETER (Size, totnnbs   );
        CONSTRUCT_PARAMETER (Size, nrays     );
        CONSTRUCT_PARAMETER (Size, hnrays    );
        CONSTRUCT_PARAMETER (Size, order_min );
        CONSTRUCT_PARAMETER (Size, order_max );
        CONSTRUCT_PARAMETER (Size, nboundary );
        CONSTRUCT_PARAMETER (Size, nfreqs    );
        CONSTRUCT_PARAMETER (Size, nspecs    );
        CONSTRUCT_PARAMETER (Size, nlspecs   );
        CONSTRUCT_PARAMETER (Size, nlines    );
        CONSTRUCT_PARAMETER (Size, nquads    );

        CONSTRUCT_PARAMETER (Real, pop_prec);

        CONSTRUCT_PARAMETER (bool, use_scattering      );
        CONSTRUCT_PARAMETER (bool, store_intensities   );
        CONSTRUCT_PARAMETER (bool, use_Ng_acceleration );
        CONSTRUCT_PARAMETER (bool, spherical_symmetry  );
        CONSTRUCT_PARAMETER (bool, adaptive_ray_tracing);
    }

    Parameters (const Parameters& parameters)
    {
        COPY_PARAMETER (model_name);

        COPY_PARAMETER (dimension );
        COPY_PARAMETER (npoints   );
        COPY_PARAMETER (totnnbs   );
        COPY_PARAMETER (nrays     );
        COPY_PARAMETER (hnrays    );
        COPY_PARAMETER (order_min );
        COPY_PARAMETER (order_max );
        COPY_PARAMETER (nboundary );
        COPY_PARAMETER (nfreqs    );
        COPY_PARAMETER (nspecs    );
        COPY_PARAMETER (nlspecs   );
        COPY_PARAMETER (nlines    );
        COPY_PARAMETER (nquads    );

        COPY_PARAMETER (pop_prec);

        COPY_PARAMETER (use_scattering      );
        COPY_PARAMETER (store_intensities   );
        COPY_PARAMETER (use_Ng_acceleration );
        COPY_PARAMETER (spherical_symmetry  );
        COPY_PARAMETER (adaptive_ray_tracing);
    }
};
