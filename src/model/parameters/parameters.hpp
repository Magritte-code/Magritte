#pragma once


#include <limits>

#include "../../configure.hpp"
#include "io/io.hpp"
#include "tools/setOnce.hpp"


#define CREATE_PARAMETER(type, x)                                               \
    private:                                                                    \
        SetOnce<type> x##__{#x};                                                \
    public:                                                                     \
        inline void set_##x (const type value)       {x##__.set(value);};\
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


    Size n_tabulated_profile_funs    = std::pow(2,14)+1;
    Size n_off_diag                  = 0;
    // diff can by computed by comparing trapezoidal rule (w*(exp-(x)^2+exp-(x+w)^2)/2) and cdf of normal (std dev=1/sqrt pi) from x to x+w
    // Real max_width_fraction          = 0.5;//max diff between the normal cdf and trapezoidal rule is +-5%
    // Real max_width_fraction          = 0.4;//max diff is +-3%
    Real max_width_fraction          = 0.35;//max diff is +-2.5%
    // Real max_width_fraction          = 0.3;//max diff is +-2%
    Real convergence_fraction        = 0.995;
    Real min_rel_pop_for_convergence = 1.0e-10;
    Real pop_prec                    = 1.0e-6;
    Real min_opacity                 = 1.0e-26;
    Real min_dtau                    = 1.0e-15;
    bool store_intensities           = false;
    bool use_Ng_acceleration         = true;
    bool one_line_approximation      = false;

    //When summing over the nearby lines, this parameter controls in what range we consider the lines to lie close
    Real max_distance_opacity_contribution = 10.0;//maximal distance at which a line can contribute to a frequency


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
