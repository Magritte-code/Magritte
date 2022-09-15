#include "paracabs.hpp"
#include "tools/types.hpp"


///  Indexer for point and line indices
///    @param[in] p          : index of the point
///    @param[in] line_index : index of the line
/////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size line_index) const
{
    return line_index + p*parameters->nlines();
}


///  Indexer for cell, line producing species and transition indices
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
////////////////////////////////////////////////////////////////////
inline Size Lines :: line_index (const Size l, const Size k) const
{
    return k + nrad_cum[l];
}


///  Indexer for point, line producing species and transition indices
///    @param[in] p : index of the point
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
/////////////////////////////////////// /////////////////////////////
inline Size Lines :: index (const Size p, const Size l, const Size k) const
{
    return index (p, line_index (l, k));
}


///  Setter for line emissivity and opacity
///////////////////////////////////////////
inline void Lines :: set_emissivity_and_opacity ()
{
    threaded_for (p, parameters->npoints(),
    {
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size lid = line_index (l, k);

                emissivity (p, lid) = lineProducingSpecies[l].get_emissivity (p, k);
                   opacity (p, lid) = lineProducingSpecies[l].get_opacity    (p, k);
            }
        }
    })
}


///  Setter for line widths
///    @param[in] thermodynamics : reference to thermodynamics module
/////////////////////////////////////////////////////////////////////
inline void Lines :: set_inverse_width (const Thermodynamics& thermodynamics)
{
    threaded_for (p, parameters->npoints(),
    {
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real invr_mass = lineProducingSpecies[l].linedata.inverse_mass;
                const Real frequency = lineProducingSpecies[l].linedata.frequency[k];

                const Size lid = line_index (l, k);

                inverse_width (p, lid) = (Real) 1.0 / thermodynamics.profile_width (invr_mass, p, frequency);
            }
        }
    })
}



/// Computes the cooling rates using the collisional approach
///   @param[in] model: the model for which to compute the cooling rates
inline void Lines :: compute_cooling_collisional(Cooling& cooling, const Double2& abundance, const Vector<Real>& temperature)
{
    // Temperature& temperature=model.thermodynamics.temperature;
    // Double2 abundance=model.chemistry.species.abundance;

    //for every point, compute the cooling rate independently
    threaded_for(p, parameters->npoints(),
    {
        Real temp_cooling_rate=0.0;
        const Real tmp = temperature[p];

        for (LineProducingSpecies &lspec : lineProducingSpecies)
        {
            //cooling is done by colliding with other species (in this formulation)
            //Implementation note: this formulation of the cooling rates is found in Sahai 1990
            //some simple setup, ensuring the collisional rates have been computed correctly
            for (CollisionPartner &colpar : lspec.linedata.colpar)
            {
                Real abn = abundance[p][colpar.num_col_partner];

                colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
                colpar.interpolate_collision_coefficients (tmp);

                for (Size k = 0; k < colpar.ncol; k++)
                {
                    const Real collision_rate_ij = colpar.Cd_intpld[k] * abn; //deexcitation->heating
                    const Real collision_rate_ji = colpar.Ce_intpld[k] * abn; //excitation  ->cooling

                    const Size I = lspec.index (p, colpar.icol[k]);
                    const Size J = lspec.index (p, colpar.jcol[k]);

                    const Real energy_diff=lspec.linedata.energy[colpar.icol[k]]-lspec.linedata.energy[colpar.jcol[k]];
                    // std::cout<<"energy_diff: "<<energy_diff<<std::endl;

                    temp_cooling_rate-=lspec.population[I]*collision_rate_ij*energy_diff;
                    temp_cooling_rate+=lspec.population[J]*collision_rate_ji*energy_diff;
                }
            }
        }
        cooling.cooling_rate[p]=temp_cooling_rate;
    });
}
