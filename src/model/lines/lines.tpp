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

///  Sets the vector (per species) of vectors (per line) of all level populations
///    @param[in]  new_population: the new level populations
//////////////////////////////////////////////////////////////////////////////////
inline void Lines :: set_all_level_pops(vector<VectorXr> new_population)
{
  threaded_for (i, parameters.nlspecs(),
    lineProducingSpecies[i].set_all_level_pops(new_population[i]);
  )
}

///  Returns the vector (per species) of vectors (per line) of all level populations
////////////////////////////////////////////////////////////////////////////////////
inline vector<VectorXr> Lines :: get_all_level_pops()
{
  vector<VectorXr> toreturn;
  toreturn.reserve(parameters.nlspecs());
  //No threaded for possible, because we might push back in a random order
  for (Size i=0; i<parameters.nlspecs(); i++)
  {
    toreturn.push_back(lineProducingSpecies[i].get_all_level_pops());
  }
  return toreturn;
}

///  Writes the level populations of a certain iteration
////////////////////////////////////////////////////////
//currently also writes the J_lin and J_eff due to using the already implemented lineProducingSpecies::write_populations
inline void Lines::write_populations_of_iteration(const Io& io, const Size it, const Size lvl) const
{
  const string tag_it="lvl"+std::to_string(lvl)+"it"+std::to_string(it);
  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      lineProducingSpecies[l].write_populations (io, l, tag_it);
  }
}

///  Reads the level populations of a given iteration
///    @param[in]:  io: Reference to the Io structure
///    @param[in]:  it: The number of the iteration to read from
///    @param[in]:  lvl: The coarsening level to read from
///////////////////////////////////////////////////////
//currently also reads the J_lin and J_eff due to using the already implemented lineProducingSpecies::write_populations
inline void Lines::read_populations_of_iteration(const Io& io, const Size it, const Size lvl)
{
  const string tag_it="lvl"+std::to_string(lvl)+"it"+std::to_string(it);
  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      lineProducingSpecies[l].read_populations (io, l, tag_it);
  }
}
