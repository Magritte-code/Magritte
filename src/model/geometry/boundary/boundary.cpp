#include "boundary.hpp"


const string prefix = "geometry/boundary/";


void Boundary :: read (const Io& io)
{
    cout << "Reading boundary..." << endl;

    // Resize boundary
    point2boundary = (size_t*) pc::accelerator::malloc (npoints*sizeof(size_t));
    is_on_boundary = (bool*)   pc::accelerator::malloc (npoints*sizeof(bool));

    // Initialise
    for (size_t p = 0; p < npoints; p++)
    {
        point2boundary[p] = npoints;
        is_on_boundary[p] = false;
    }

    // Read boundary list
    io.read_length (prefix+"boundary2point", nboundary);

    Long1 boundary2point_buffer (nboundary);

    boundary2point = (size_t*) pc::accelerator::malloc (nboundary*sizeof(size_t));

    io.read_list (prefix+"boundary2point", boundary2point_buffer);

    // Set helper variables to identify the boundary
    for (size_t b = 0; b < nboundary; b++)
    {
        boundary2point[b]                 = boundary2point_buffer[b];
        point2boundary[boundary2point[b]] = b;
        is_on_boundary[boundary2point[b]] = true;
    }

    cout << "n boundary = " << nboundary << endl;
}


void Boundary :: write (const Io& io) const
{
    Long1 boundary2point_buffer (nboundary);

    for (size_t b = 0; b < nboundary; b++)
    {
        boundary2point_buffer[b] = boundary2point[b];
    }

    io.write_list (prefix+"boundary2point", boundary2point_buffer);
}