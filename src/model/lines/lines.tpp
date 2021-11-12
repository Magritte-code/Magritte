#include "paracabs.hpp"
#include "tools/types.hpp"


///  Indexer for point and line indices
///    @param[in] p          : index of the point
///    @param[in] line_index : index of the line
/////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size line_index) const
{
    return line_index + p*parameters.nlines();
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
    threaded_for (p, parameters.npoints(),
    {
        for (Size l = 0; l < parameters.nlspecs(); l++)
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
    threaded_for (p, parameters.npoints(),
    {
        for (Size l = 0; l < parameters.nlspecs(); l++)
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


///  Calculates the cooling rates
///  NOTE: call this only after computing the radiation field
///////////////////////////////////////////
inline void Lines :: calculate_cooling_rates (const Double2 &abundance, const Vector<Real> &temperature)
{
    // Matrix<Real> partial_cooling_rates; //for each point and line species
    // partial_cooling_rates.resize (parameters.npoints(),  parameters.nlspecs());
    //FIXME: change the colpar.adjust_abundance... noncense to finally allow parallelisation ;-; same for interpolate_collision_coefficients
    // threaded_for (p, parameters.npoints(),
    for (Size p=0; p<parameters.npoints(); p++)
    {
        Real cooling_of_point=0;

        const Real tmp = temperature[p];

        // TODO: if certain species may not participate, exclude them here
        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            Real cooling_of_lspec=0;

            for (CollisionPartner &colpar : lineProducingSpecies[l].linedata.colpar)
            {
                Real cooling_with_colpar=0;

                Real abn = abundance[p][colpar.num_col_partner];

                colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
                colpar.interpolate_collision_coefficients (tmp);

                for (Size k = 0; k < colpar.ncol; k++)
                {
                    const Real collision_rate_ij = colpar.Cd_intpld[k] * abn; //deexcitation->heating
                    const Real collision_rate_ji = colpar.Ce_intpld[k] * abn; //excitation  ->cooling

                    const Size I = lineProducingSpecies[l].index (p, colpar.icol[k]);
                    const Size J = lineProducingSpecies[l].index (p, colpar.jcol[k]);

                    const Real energy_diff=lineProducingSpecies[l].linedata.energy[colpar.icol[k]]-lineProducingSpecies[l].linedata.energy[colpar.jcol[k]];
                    // std::cout<<"energy_diff: "<<energy_diff<<std::endl;

                    cooling_with_colpar-=lineProducingSpecies[l].population[I]*collision_rate_ij*energy_diff;
                    cooling_with_colpar+=lineProducingSpecies[l].population[J]*collision_rate_ji*energy_diff;

                }

                cooling_of_lspec+=cooling_with_colpar;
            }

            // for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            // {
            //     const Size lid = line_index (l, k);
            //
            //     //TODO: replace this with deltaE(=hv)*Collisional rates()(*angle?)
            //     //This calculates 'cooling' by computing what energy escapes by radiation (not by ckecking what energy is extracted from the gas (towards quantum states))
            //     //This is can be wrong and should be replaced
            //     cooling_of_lspec+=emissivity(p, lid)*line[lid]*FOUR_PI;//emissivity*freq(->energy)*4pi(->integrate out direction)
            //     cooling_of_lspec-=opacity   (p, lid)*line[lid]*FOUR_PI*lineProducingSpecies[l].Jlin[p][k];//also times averaged intensity
            // }

            cooling_of_point+=cooling_of_lspec;
            //TODO: figure out whether we need the individual cooling rates of the line producing species
        }
        // //TODO: add fancy reduction operation instead of this
        cooling_rates[p]=cooling_of_point;
    }//)

}
