#pragma once


#include <limits>

#include "../../configure.hpp"
#include "io/io.hpp"
#include "tools/setOnce.hpp"


#define CREATE_PARAMETER(type, x)                                               \
    private:                                                                    \
        SetOnce<type> x##__;                                                    \
    public:                                                                     \
        inline void set_##x (const type value)       {       x##__.set(value);};\
        inline type       x (                ) const {return x##__.get(     );};


///  Parameters: secure structure for the model parameters
//////////////////////////////////////////////////////////
struct Parameters
{
    CREATE_PARAMETER (string, version);

    CREATE_PARAMETER (string, model_name);

    CREATE_PARAMETER (Size, dimension );
    CREATE_PARAMETER (Size, npoints   );
    CREATE_PARAMETER (Size, nrays     );
    CREATE_PARAMETER (Size, hnrays    );
    CREATE_PARAMETER (Size, nboundary );
    CREATE_PARAMETER (Size, nfreqs    );
    CREATE_PARAMETER (Size, nspecs    );
    CREATE_PARAMETER (Size, nlspecs   );
    CREATE_PARAMETER (Size, nlines    );
    CREATE_PARAMETER (Size, nquads    );

    CREATE_PARAMETER (bool, use_scattering        );
    CREATE_PARAMETER (bool, spherical_symmetry    );
    CREATE_PARAMETER (bool, adaptive_ray_tracing  );


    Size n_off_diag                  = 0;
    Real max_width_fraction          = 0.5;
    Real convergence_fraction        = 0.995;
    Real min_rel_pop_for_convergence = 1.0e-10;
    Real pop_prec                    = 1.0e-6;
    Real min_opacity                 = 1.0e-26;
    bool store_intensities           = false;
    bool use_Ng_acceleration         = true;
    bool one_line_approximation      = false;


    Parameters ()
    {
        // Disable scattering
        use_scattering__.set (false);

        // Set version
        version__.set (MAGRITTE_VERSION);

        // Set defaults
        spherical_symmetry__  .set_default (false);
        adaptive_ray_tracing__.set_default (false);
    };


    void read (const Io &io);
    void write(const Io &io) const;
};
