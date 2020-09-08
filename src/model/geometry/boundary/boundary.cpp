#include "boundary.hpp"


const string prefix = "geometry/boundary/";


void Boundary :: read (const Io& io)
{
    cout << "Reading boundary..." << endl;

    // Resize boundary
    point2boundary.resize (parameters.npoints());

    // Initialise
    for (Size p = 0; p < parameters.npoints(); p++)
    {
        point2boundary[p] = parameters.npoints();
    }

    // Read boundary list
    parameters.set_nboundary (io.get_length (prefix+"boundary2point"));

    boundary2point.resize (parameters.nboundary());

    io.read_list (prefix+"boundary2point", boundary2point);

    // Set helper variables to identify the boundary
    for (Size b = 0; b < parameters.nboundary(); b++)
    {
        point2boundary[boundary2point[b]] = b;
    }

    cout << "n boundary = " << parameters.nboundary() << endl;

    boundary2point.copy_vec_to_ptr ();
    point2boundary.copy_vec_to_ptr ();
}


void Boundary :: write (const Io& io) const
{
    io.write_list (prefix+"boundary2point", boundary2point);
}