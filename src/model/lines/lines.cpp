#include "lines.hpp"
#include "tools/heapsort.hpp"


const string prefix = "lines/";


///  Reader for the Lines data
///    @param[in] io         : io object to read with
///    @param[in] parameters : model parameters object
//////////////////////////////////////////////////////
void Lines :: read (const Io& io)
{
    cout << "Reading lines..." << endl;

    /// Read and set nlspecs
    parameters.set_nlspecs (io.get_length (prefix+"lineProducingSpecies_"));

    /// Read line producing species data
    lineProducingSpecies.resize (parameters.nlspecs());

    for (Size l = 0; l < parameters.nlspecs(); l++)
    {
        lineProducingSpecies[l].read (io, l);
    }

    /// Set nrad_cum, a helper variable for determining indices
    nrad_cum.resize (parameters.nlspecs());
    nrad_cum[0] = 0;

    for (Size l = 1; l < parameters.nlspecs(); l++)
    {
        nrad_cum[l] = nrad_cum[l-1] + lineProducingSpecies[l-1].linedata.nrad;
    }

    /// Extract nlines
    Size nlines = 0;

    for (const LineProducingSpecies &lspec : lineProducingSpecies)
    {
        nlines += lspec.linedata.nrad;
    }

    parameters.set_nlines (nlines);

    /// Set and sort lines and their indices
    line      .resize (parameters.nlines());
    // line_index.resize (parameters.nlines());

    Size index = 0;

    for (const LineProducingSpecies &lspec : lineProducingSpecies)
    {
        for (Size k = 0; k < lspec.linedata.nrad; k++)
        {
            line      [index] = lspec.linedata.frequency[k];
            // line_index[index] = index;
            index++;
        }
    }


    // heapsort (line, line_index);


    // emissivity.resize (parameters.npoints()*parameters.nlines());
    // opacity   .resize (parameters.npoints()*parameters.nlines());

    emissivity   .resize (parameters.npoints(), parameters.nlines());
    opacity      .resize (parameters.npoints(), parameters.nlines());
    inverse_width.resize (parameters.npoints(), parameters.nlines());
}


///  Writer for the Lines data
///    @param[in] io: io object to write with
/////////////////////////////////////////////
void Lines :: write (const Io& io) const
{
    cout << "Writing lines..." << endl;

    for (Size l = 0; l < lineProducingSpecies.size(); l++)
    {
        lineProducingSpecies[l].write (io, l);
    }
}


void Lines :: iteration_using_LTE (const Double2 &abundance, const Vector<Real> &temperature)
{
    for (LineProducingSpecies &lspec : lineProducingSpecies)
    {
        lspec.update_using_LTE (abundance, temperature);
    }

    set_emissivity_and_opacity ();

    //gather_emissivities_and_opacities ();
}


void Lines :: iteration_using_Ng_acceleration (const Real pop_prec, vector<Size> &points_in_grid)
{
    for (LineProducingSpecies &lspec : lineProducingSpecies)
    {
        lspec.update_using_Ng_acceleration ();
        lspec.check_for_convergence        (pop_prec,points_in_grid);
    }

    set_emissivity_and_opacity ();

    //gather_emissivities_and_opacities ();
}


void Lines :: iteration_using_statistical_equilibrium (
    const Double2      &abundance,
    const Vector<Real> &temperature,
    const Real          pop_prec,
    vector<Size> &points_in_grid)
{
    for (LineProducingSpecies &lspec : lineProducingSpecies)
    {
        lspec.update_using_statistical_equilibrium (abundance, temperature, points_in_grid);
        std::cout<<"checking for convergence"<<std::endl;
        lspec.check_for_convergence                (pop_prec, points_in_grid);
    }
    std::cout<<"setting emissivity and opacity"<<std::endl;
    set_emissivity_and_opacity ();

    //gather_emissivities_and_opacities ();
}


void Lines :: gather_emissivities_and_opacities ()

#if (MPI_PARALLEL)

{
    // Get number of processes
    const int comm_size = MPI_comm_size ();

    // Extract the buffer lengths and displacements
    int *buffer_lengths = new int[comm_size];
    int *displacements  = new int[comm_size];

    for (int w = 0; w < comm_size; w++)
    {
        long start = ( w   *parameters.npoints())/comm_size;
        long stop  = ((w+1)*parameters.npoints())/comm_size;

        long ncells_red_w = stop - start;

        buffer_lengths[w] = ncells_red_w * parameters.nlines();
    }

    displacements[0] = 0;

    for (int w = 1; w < comm_size; w++)
    {
        displacements[w] = buffer_lengths[w-1];
    }

    // Call MPI to gather the emissivity data
    int ierr_em = MPI_Allgatherv (
                      MPI_IN_PLACE,            // pointer to data to be send (here in place)
                      0,                       // number of elements in the send buffer
                      MPI_DATATYPE_NULL,       // type of the send data
                      emissivity.data(),       // pointer to the data to be received
                      buffer_lengths,          // number of elements in receive buffer
                      displacements,           // displacements between data blocks
  	                  MPI_DOUBLE,              // type of the received data
                      MPI_COMM_WORLD );
    assert (ierr_em == 0);

    // Call MPI to gather the opacity data
    int ierr_op = MPI_Allgatherv (
                      MPI_IN_PLACE,            // pointer to data to be send (here in place)
                      0,                       // number of elements in the send buffer
                      MPI_DATATYPE_NULL,       // type of the send data
                	  opacity.data(),          // pointer to the data to be received
                	  buffer_lengths,          // number of elements in receive buffer
                      displacements,           // displacements between data blocks
                	  MPI_DOUBLE,              // type of the received data
                      MPI_COMM_WORLD );
    assert (ierr_op == 0);

    delete [] buffer_lengths;
    delete [] displacements;
}

#else

{
    return;
}

#endif
