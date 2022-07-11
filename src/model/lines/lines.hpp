#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Lines
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    vector <LineProducingSpecies> lineProducingSpecies;

    Vector<Real> line;                       ///< [Hz] line center frequencies (NOT ordered!)
    Vector<Size> nrad_cum;                   ///< Cumulative number of radiative transitions

    Matrix<Real> emissivity;                 ///< line emissivity    (p, lid)
    Matrix<Real> opacity;                    ///< line opacity       (p, lid)
    Matrix<Real> inverse_width;              ///< inverse line width (p, lid)


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

    void resize_LineProducingSpecies(const Size nlspec)
    {
        lineProducingSpecies.resize (parameters->nlspecs(), LineProducingSpecies(parameters));
    }

    void PORTAL_solve_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature );

    void gather_emissivities_and_opacities ();
};


#include "lines.tpp"
