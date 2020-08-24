#include "lineProducingSpecies.hpp"


const string prefix = "lines/lineProducingSpecies_";


///  Reader for the LineProducingSpecies data from the Io object
///    @param[in] io : io object
///    @param[in] l  : nr of line producing species
////////////////////////////////////////////////////////////////
void LineProducingSpecies :: read (const Io& io, const Size l)
{
    cout << "Reading lineProducingSpecies..." << endl;

    linedata  .read (io, l);
    quadrature.read (io, l);

    nquads = quadrture.get_nquads();

    RT        .resize (npoints*linedata.nlev, npoints*linedata.nlev);
    LambdaStar.resize (npoints*linedata.nlev, npoints*linedata.nlev);
    LambdaTest.resize (npoints*linedata.nlev, npoints*linedata.nlev);

    lambda.initialize (npoints, linedata.nrad);

    Jeff.resize (npoints);
    Jlin.resize (npoints);
    Jdif.resize (npoints);

    for (Size p = 0; p < npoints; p++)
    {
        Jeff[p].resize (linedata.nrad);
        Jlin[p].resize (linedata.nrad);
        Jdif[p].resize (linedata.nrad);
    }

    nr_line.resize (npoints);

    for (Size p = 0; p < npoints; p++)
    {
        nr_line[p].resize (linedata.nrad);

        for (Size k = 0; k < linedata.nrad; k++)
        {
            nr_line[p][k].resize (nquads);
        }
    }

    population_prev1.resize (npoints*linedata.nlev);
    population_prev2.resize (npoints*linedata.nlev);
    population_prev3.resize (npoints*linedata.nlev);
      population_tot.resize (npoints*linedata.nlev);
          population.resize (npoints*linedata.nlev);

    const string prefix_l = prefix + std::to_string (l) + "/";

    io.read_list (prefix_l+"population_tot", population_tot);

    read_populations (io, l, "");

    Double2 pops_prev1 (npoints, Double1 (linedata.nlev));
    Double2 pops_prev2 (npoints, Double1 (linedata.nlev));
    Double2 pops_prev3 (npoints, Double1 (linedata.nlev));

    int err_prev1 = io.read_array (prefix_l+"population_prev1", pops_prev1);
    int err_prev2 = io.read_array (prefix_l+"population_prev2", pops_prev2);
    int err_prev3 = io.read_array (prefix_l+"population_prev3", pops_prev3);

    threaded_for (p, ncells)
    {
        for (Size i = 0; i < linedata.nlev; i++)
        {
            if (err_prev1 == 0) {population_prev1 (index (p, i)) = pops_prev1[p][i];}
            if (err_prev2 == 0) {population_prev2 (index (p, i)) = pops_prev2[p][i];}
            if (err_prev3 == 0) {population_prev3 (index (p, i)) = pops_prev3[p][i];}
        }
    }
}




///  Writer for the LineProducingSpecies data to the Io object
///    @param[in] io : io object
///    @param[in] l  : nr of line producing species
//////////////////////////////////////////////////////////////
void LineProducingSpecies :: write (const Io &io, const long l) const
{
    cout << "Writing lineProducingSpecies..." << endl;

    linedata.write (io, l);

    quadrature.write (io, l);

    write_populations (io, l, "");

    const string prefix_l = prefix + std::to_string (l) + "/";

    io.write_list (prefix_l+"population_tot", population_tot);

    Double2 pops_prev1 (npoints, Double1 (linedata.nlev));
    Double2 pops_prev2 (npoints, Double1 (linedata.nlev));
    Double2 pops_prev3 (npoints, Double1 (linedata.nlev));

    threaded_for (p, npoints)
    {
        for (Size i = 0; i < linedata.nlev; i++)
        {
            pops_prev1[p][i] = population_prev1 (index (p, i));
            pops_prev2[p][i] = population_prev2 (index (p, i));
            pops_prev3[p][i] = population_prev3 (index (p, i));
        }
    }

    io.write_array (prefix_l+"population_prev1", pops_prev1);
    io.write_array (prefix_l+"population_prev2", pops_prev2);
    io.write_array (prefix_l+"population_prev3", pops_prev3);





///  Reader for the level populations from the Io object
///    @param[in] io  : io object
///    @param[in] l   : number of line producing species
///    @param[in] tag : extra info tag
////////////////////////////////////////////////////////

int LineProducingSpecies :: read_populations (const Io &io, const long l, const string tag)
{
    const string prefix_l = prefix + std::to_string (l) + "/";

    Double2 pops (ncells, Double1 (linedata.nlev));

    int err = io.read_array (prefix_l+"population"+tag, pops);

    if (err == 0)
    {
        OMP_PARALLEL_FOR (p, ncells)
        {
            for (long i = 0; i < linedata.nlev; i++)
            {
                population (index (p, i)) = pops[p][i];
            }
        }
    }

    io.read_array (prefix_l+"J_lin"+tag, Jlin);
    io.read_array (prefix_l+"J_eff"+tag, Jeff);

    return (0);
}




///  Writer for the level populations to the Io object
///    @param[in] io  : io object
///    @param[in] l   : number of line producing species
///    @param[in] tag : extra info tag
////////////////////////////////////////////////////////

int LineProducingSpecies :: write_populations (const Io &io, const long l, const string tag) const
{
  const string prefix_l = prefix + std::to_string (l) + "/";

  Double2 pops (ncells, Double1 (linedata.nlev));

  cout << "Writing populations..." << endl;

  OMP_PARALLEL_FOR (p, ncells)
  {
      for (long i = 0; i < linedata.nlev; i++)
      {
          pops[p][i] = population (index (p, i));
      }
  }

  io.write_array (prefix_l+"population"+tag, pops);

  io.write_array (prefix_l+"J_lin"+tag, Jlin);
  io.write_array (prefix_l+"J_eff"+tag, Jeff);

  return (0);
}
