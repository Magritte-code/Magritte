#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Lines
{
    Parameters parameters;

    vector <LineProducingSpecies> lineProducingSpecies;

    // Real1 line;         ///< [Hz] line center frequencies orderd
    Vector<Real> line;     ///< [Hz] line center frequencies (NOT ordered!)
    // Size1 line_index;   ///< index of the corresponding frequency in line

    Vector<Size> nrad_cum;

    // Real1 emissivity;   ///< line emissivity (p,l,k)
    // Real1 opacity;      ///< line opacity    (p,l,k)


    // Cooling rates can also be moved somewhere else, but seems logical here
    Vector<Real> cooling_rates;   ///<the computed cooling rates

    Matrix<Real> emissivity;      ///< line emissivity    (p, lid)
    Matrix<Real> opacity;         ///< line opacity       (p, lid)
    Matrix<Real> inverse_width;   ///< inverse line width (p, lid)

    void read  (const Io& io);
    void write (const Io& io) const;

    void iteration_using_LTE (
        const Double2      &abundance,
        const Vector<Real> &temperature);

    void iteration_using_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature,
        const Real          pop_prec             );

    void iteration_using_Ng_acceleration (
        const Real pop_prec              );

    inline Size      index (const Size p, const Size line_index     ) const;
    inline Size line_index (              const Size l, const Size k) const;
    inline Size      index (const Size p, const Size l, const Size k) const;

    inline void set_emissivity_and_opacity ();
    inline void set_inverse_width (const Thermodynamics& thermodynamics);

    void gather_emissivities_and_opacities ();

    // Calculates the cooling rates for all points, currently for all species
    inline void calculate_cooling_rates ();
};


#include "lines.tpp"
