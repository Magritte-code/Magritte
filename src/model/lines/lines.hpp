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

    Vector<Real> line;                     ///< [Hz] line center frequencies (NOT ordered!)
    Matrix<Real> line_frequency;           ///< line frequency
    Vector<Real> line_inverse_mass;        ///< inverse mass 
    Matrix<Real> line_A;                   ///< Einstein A coefficient
    Matrix<Real> line_quadrature_weight;   ///< quadrature weights

    Vector<Size> nrad_cum;

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
};


#include "lines.tpp"
