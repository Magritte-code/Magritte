#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"
#include "model/thermodynamics/thermodynamics.hpp"
#include "model/cooling/cooling.hpp"
#include "model/radiation/radiation.hpp"


struct Lines
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    vector <LineProducingSpecies> lineProducingSpecies;

    Vector<Real> line;                       ///< [Hz] line center frequencies (NOT ordered!)
    Vector<Size> nrad_cum;                   ///< Cumulative number of radiative transitions

    Real1 sorted_line;                       ///< [Hz] sorted line center frequencies (ordered)
    Size1 sorted_line_map;                   ///< Mapping from sorted line indices to the default ones
    Real max_inverse_mass;                   ///< [1/AMU] maximal inverse mass of all line species

    Matrix<Real> emissivity;                 ///< line emissivity    (p, lid)
    Matrix<Real> opacity;                    ///< line opacity       (p, lid)
    Matrix<Real> inverse_width;              ///< inverse line width (p, lid)

    Matrix<Real> line_cooling_rate;          ///< cooling rate per line (p, lid) UNITS TODO
    Vector<Real> total_cooling_rate;         ///< cooling rate per point (p)


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

    //Compute the cooling rates using the collisional formulation
    void compute_cooling_collisional(Cooling& cooling, const Double2& abundance, const Vector<Real>& temperature);
    void compute_line_cooling_radiative(Cooling& cooling, Radiation& radiation);
    void compute_line_cooling_radiative_sparse(Cooling& cooling);

    void gather_emissivities_and_opacities ();
};


#include "lines.tpp"
