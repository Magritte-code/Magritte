#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Lines
{
    const Size LOG2_N_TABULED_PROFILE_FUNS_M1=14;//=2^14+1
    const Size N_TABULATED_PROFILE_FUNS = std::pow(2, LOG2_N_TABULED_PROFILE_FUNS_M1) + 1;
    const Size LOG2_MAX_DISTANCE_INTERVAL=4;//=2^4=16
    const double MAX_DISTANCE_INTERVAL = std::pow(2, LOG2_MAX_DISTANCE_INTERVAL-1);//=2^3
    const Size NET_LEFT_SHIFT = LOG2_N_TABULED_PROFILE_FUNS_M1 - LOG2_MAX_DISTANCE_INTERVAL;
    const Size MULTIPLICATION_FACTOR = std::pow(2, NET_LEFT_SHIFT);

    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    vector <LineProducingSpecies> lineProducingSpecies;

    Vector<Real> line;                       ///< [Hz] line center frequencies (NOT ordered!)
    Vector<Size> nrad_cum;                   ///< Cumulative number of radiative transitions

    Matrix<Real> emissivity;                 ///< line emissivity    (p, lid)
    Matrix<Real> opacity;                    ///< line opacity       (p, lid)
    Matrix<Real> inverse_width;              ///< inverse line width (p, lid)


    Vector<Real> tabulated_gaussians;         ///< Table containing tabulated values for gaussians (the e^(-x^2) part)
                                              /// Computed at each x in a linear interval [-parameters.max_distance_opacity_contribution, +parameters.max_distance_opacity_contribution]
    Vector<Real> tabulated_erfs;              ///< Table containing tabulated values for gaussians (the e^(-x^2) part)
                                              /// Computed at each x in a linear interval [-parameters.max_distance_opacity_contribution, +parameters.max_distance_opacity_contribution]


    Lines (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;

    void iteration_using_LTE (
        const Double2      &abundance,
        const Vector<Real> &temperature);

    void iteration_using_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature,
        const Real          pop_prec             );

    void iteration_using_statistical_equilibrium_sparse (
        const Double2      &abundance,
        const Vector<Real> &temperature,
        const Real          pop_prec                    );

    void iteration_using_Ng_acceleration (
        const Real pop_prec              );

    inline Size      index (const Size p, const Size line_index     ) const;
    inline Size line_index (              const Size l, const Size k) const;
    inline Size      index (const Size p, const Size l, const Size k) const;

    inline void set_emissivity_and_opacity ();
    inline void set_inverse_width (const Thermodynamics& thermodynamics);

    inline Size convert_to_table_index (double x) const;
    inline void set_tabulated_gaussians ();
    inline void set_tabulated_erfs ();
    inline Real compute_tabulated_gaussian (double x) const;
    inline Real compute_tabulated_erf (double x) const;

    void resize_LineProducingSpecies(const Size nlspec)
    {
        lineProducingSpecies.resize (parameters->nlspecs(), LineProducingSpecies(parameters));
    }

    void gather_emissivities_and_opacities ();
};


#include "lines.tpp"
