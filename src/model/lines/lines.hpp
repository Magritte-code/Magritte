#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"


struct Lines
{
    Parameters parameters;

    std::vector <LineProducingSpecies> lineProducingSpecies;

    Real1 line;         ///< [Hz] line center frequencies orderd
    Size1 line_index;   ///< index of the corresponding frequency in line

    Size1 nrad_cum;

    Real1 emissivity;   ///< line emissivity (p,l,k)
    Real1 opacity;      ///< line opacity    (p,l,k)

    void read  (const Io& io);
    void write (const Io& io) const;

    void iteration_using_LTE (
        const Double2 &abundance,
        const Real1   &temperature);

    void iteration_using_statistical_equilibrium (
        const Double2 &abundance,
        const Real1   &temperature,
        const Real     pop_prec                  );

    void iteration_using_Ng_acceleration (
        const Real pop_prec              );


    inline Size index (const Size p, const Size l, const Size k) const;
    inline Size index (const Size p, const Size line_index     ) const;

    inline void set_emissivity_and_opacity ();

    void gather_emissivities_and_opacities ();
};


#include "lines.tpp"