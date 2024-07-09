#pragma once

#include "../../configure.hpp"
#include "io/io.hpp"
#include "tools/setOnce.hpp"

#include <limits>

#define CREATE_PARAMETER(type, x, essential)                     \
  private:                                                       \
    SetOnce<type, essential> x##__{#x};                          \
                                                                 \
  public:                                                        \
    inline void set_##x(const type value) { x##__.set(value); }; \
    inline type x() const { return x##__.get(); };

// Global variables may change their value at fixed times,
// explained in the comments, in constrast to static
// parameters
#define CREATE_MUTABLE_PARAMETER(type, x)                     \
  private:                                                    \
    type x##__;                                               \
                                                              \
  public:                                                     \
    inline void set_##x(const type value) { x##__ = value; }; \
    inline type x() const { return x##__; };

// template<typename T>
// CREATE_GLOBAL_VARIABLE(string name)
// {
//     private:
//     T x
// }

///  Parameters: secure structure for the model parameters
//////////////////////////////////////////////////////////
struct Parameters {
    // essential parameters which need to be set by the user
    CREATE_PARAMETER(string, version, true);
    CREATE_PARAMETER(string, model_name, true);

    CREATE_PARAMETER(Size, dimension, true);
    CREATE_PARAMETER(Size, npoints, true);
    CREATE_PARAMETER(Size, nrays, true);
    CREATE_PARAMETER(Size, nboundary, true);
    CREATE_PARAMETER(Size, nspecs, true);
    CREATE_PARAMETER(Size, nlspecs, true);
    CREATE_PARAMETER(Size, nquads, true);

    CREATE_PARAMETER(bool, spherical_symmetry, true);

    // nonessential parameters which are automatically
    // inferred by magritte
    CREATE_PARAMETER(Size, hnrays, false);
    CREATE_PARAMETER(Size, nlines, false);

    // parameters which can be redefined
    CREATE_MUTABLE_PARAMETER(Size, nfreqs); // due to the fact that line
                                            // radiative transfer and the imager
                                            // require different amounts of

    // Old things that are no more connected to anything...
    CREATE_PARAMETER(bool, use_scattering, false);
    CREATE_PARAMETER(bool, adaptive_ray_tracing, false);

    Size n_off_diag = 0;
    // diff can by computed by comparing trapezoidal rule
    // (w*(exp-(x)^2+exp-(x+w)^2)/2) and cdf of normal (std
    // dev=1/sqrt pi) from x to x+w Real max_width_fraction
    // = 0.5;//max diff between the normal cdf and
    // trapezoidal rule is +-5% Real max_width_fraction =
    // 0.4;//max diff is +-3%
    Real max_width_fraction = 0.35; // max diff is +-2.5%
    // Real max_width_fraction          = 0.3;//max diff is
    // +-2%
    Real convergence_fraction        = 0.995;
    Real min_rel_pop_for_convergence = 1.0e-10;
    Real pop_prec                    = 1.0e-6;
    Real min_opacity                 = 1.0e-26;
    long double min_line_opacity     = 1.0e-13;
    Real min_dtau                    = 1.0e-15;
    Real population_inversion_fraction =
        -1.0; // previously 1.01; // threshold factor for population inversion required
              //  for LTE to be used; set this higher than 1
    bool store_intensities = false;
    // bool use_Ng_acceleration         = true;//Not used,
    // so may safely be removed
    bool use_adaptive_Ng_acceleration = true;    // whether to use an adaptive version of Ng
                                                 // acceleration; only relevant if
                                                 // use_Ng_acceleration = true
    Size Ng_acceleration_mem_limit = 9;          // determines how many previous iterations we
                                                 // hold in memory when using adaptive
                                                 // ng-acceleration
    Size adaptive_Ng_acceleration_min_order = 2; // Minimal order of adaptive ng acceleration
                                                 // used. Has to be larger than 1.
    bool adaptive_Ng_acceleration_use_max_criterion =
        true;                              // Whether or not to use max change as
                                           // criterion for adaptive ng acceleration;
                                           // uses mean change if false
    Size Ng_acceleration_remove_N_its = 0; // Number of iterations to throw away
                                           // when using adaptive ng acceleration

    /// Approximations for summing over lines; by default,
    /// we only sum over the close lines in order to compute
    /// opacity/emissivity (leading to O(Nline * ln(Nlines))
    /// scaling behavior)

    // Warning: although this has O(Nlines) scaling
    // behavior, one must check that no lines can overlap
    bool one_line_approximation = false;

    // Warning: this is has O(Nlines^2) scaling behavior, so
    // it is slow when computing using a lot of lines Note:
    // this is just legacy behavior, for those who wish to
    // compare the obtained results to previous results
    bool sum_opacity_emissivity_over_all_lines = false;

    // When summing over the nearby lines, this parameter
    // controls in what range we consider the lines to lie
    // close
    Real max_distance_opacity_contribution = 10.0; // maximal distance at which a line can
                                                   // contribute to a frequency

    Parameters() {
        // Disable scattering
        use_scattering__.set(false);

        // Set version
        version__.set(MAGRITTE_VERSION);

        // Set defaults
        spherical_symmetry__.set_default(false);
        adaptive_ray_tracing__.set_default(false);
    };

    void read(const Io& io);
    void write(const Io& io) const;
};
