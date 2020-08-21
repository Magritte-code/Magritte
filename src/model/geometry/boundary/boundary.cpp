#include "boundary.hpp"


const string prefix = "geometry/boundary/";


void Boundary :: read (const Io& io)
{
    cout << "Reading boundary..." << endl;

    // Resize boundary
    point2boundary.resize (npoints);
//    is_on_boundary.resize (npoints);

//    point2boundary = (size_t*) pc::accelerator::malloc (npoints*sizeof(size_t));
//    is_on_boundary = (bool*)   pc::accelerator::malloc (npoints*sizeof(bool));

    // Initialise
    for (Size p = 0; p < npoints; p++)
    {
        point2boundary.vec[p] = npoints;
//        is_on_boundary.vec[p] = false;
    }

    // Read boundary list
    io.read_length (prefix+"boundary2point", nboundary);


//    Long1 boundary2point_buffer (nboundary);

    boundary2point.resize (nboundary);
//    boundary2point = (size_t*) pc::accelerator::malloc (nboundary*sizeof(size_t));

    io.read_list (prefix+"boundary2point", boundary2point.vec);
//    io.read_list (prefix+"boundary2point", boundary2point_buffer);

//    nboundary++;
//    boundary2point.vec.push_back (820469);

    // Set helper variables to identify the boundary
    for (Size b = 0; b < nboundary; b++)
    {
//        boundary2point[b]                 = boundary2point_buffer[b];
        point2boundary.vec[boundary2point.vec[b]] = b;
//        is_on_boundary.vec[boundary2point.vec[b]] = true;
    }

    cout << "n boundary = " << nboundary << endl;

    cout << "copying..." << endl;
    boundary2point.copy_vec_to_ptr ();
    cout << "copying still..." << endl;

    point2boundary.copy_vec_to_ptr ();
    cout << "copying done" << endl;
//    is_on_boundary.copy_vec_to_ptr ();
}


void Boundary :: write (const Io& io) const
{
//    Long1 boundary2point_buffer (nboundary);
//
//    for (size_t b = 0; b < nboundary; b++)
//    {
//        boundary2point_buffer[b] = boundary2point[b];
//    }

//    io.write_list (prefix+"boundary2point", boundary2point_buffer);
    io.write_list (prefix+"boundary2point", boundary2point.vec);
}