#include <assert.h>
#include "points.hpp"


const string prefix = "geometry/points/";


void Points :: read (const Io& io)
{
    cout << "Reading points..." << endl;

    Size n1, n2;
    io.read_length (prefix+"position", n1);
    io.read_length (prefix+"velocity", n2);
    assert (n1 == n2);
    npoints = n1;

    cout << "npoints = " << npoints << endl;

    position.resize (npoints);
    velocity.resize (npoints);

//    position = (Vector3D*) pc::accelerator::malloc (npoints*sizeof(Vector3D));
//    velocity = (Vector3D*) pc::accelerator::malloc (npoints*sizeof(Vector3D));

    Double2 position_buffer (npoints, Double1(3));
    Double2 velocity_buffer (npoints, Double1(3));

    io.read_array (prefix+"position", position_buffer);
    io.read_array (prefix+"velocity", velocity_buffer);

    for (Size p = 0; p < npoints; p++)
    {
        position.vec[p] = Vector3D (position_buffer[p][0],
                                    position_buffer[p][1],
                                    position_buffer[p][2] );

//        position.vec[p].print();
//        cout << position.vec[p].x() << endl;

        velocity.vec[p] = Vector3D (velocity_buffer[p][0],
                                    velocity_buffer[p][1],
                                    velocity_buffer[p][2] );
    }

    io.read_length (prefix+"neighbors", tot_n_neighbors);

    cout << "tot_n_neighbors = " << tot_n_neighbors << endl;

    cum_n_neighbors.resize (npoints);
        n_neighbors.resize (npoints);
          neighbors.resize (tot_n_neighbors);

//    cum_n_neighbors = (size_t*) pc::accelerator::malloc (npoints        *sizeof(size_t));
//        n_neighbors = (size_t*) pc::accelerator::malloc (npoints        *sizeof(size_t));
//          neighbors = (size_t*) pc::accelerator::malloc (tot_n_neighbors*sizeof(size_t));

    cout << "memory_allocated = " << endl;

//    Long1 n_neighbors_buffer (npoints);
//    Long1   neighbors_buffer (tot_n_neighbors);

    cout << "lists made = " << endl;

    io.read_list (prefix+"n_neighbors", n_neighbors.vec);
    io.read_list (prefix+  "neighbors",   neighbors.vec);

    cout << "lists read = " << endl;

//        n_neighbors[0] = n_neighbors_buffer[0];
    cum_n_neighbors.vec[0] = 0;

    cout << "first put" << endl;

    for (Size p = 1; p < npoints; p++)
    {
//            n_neighbors[p] = n_neighbors_buffer[p];
        cum_n_neighbors.vec[p] = cum_n_neighbors.vec[p-1] + n_neighbors.vec[p-1];
    }

    cout << "points put" << endl;
//    cout << "n " << neighbors_buffer.size() << endl;

//    for (size_t p = 0; p < tot_n_neighbors; p++)
//    {
//        cout << p << endl;
//        neighbors[p] = neighbors_buffer[p];
//    }

    position.copy_vec_to_ptr ();
    velocity.copy_vec_to_ptr ();

    cum_n_neighbors.copy_vec_to_ptr ();
        n_neighbors.copy_vec_to_ptr ();
          neighbors.copy_vec_to_ptr ();

    cout << "neighbors put" << endl;



    nbs.resize (npoints*nnbs);

    for (Size p = 0; p < npoints; p++)
    {
        const Size     n_nbs =     n_neighbors.vec[p];
        const Size cum_n_nbs = cum_n_neighbors.vec[p];

        for (Size i = 0; (i < n_nbs) && (i < nnbs); i++)
        {
            nbs.vec[p*nnbs+i] = neighbors.vec[cum_n_nbs+i];
        }
        for (Size i = n_nbs; i < nnbs; i++)
        {
            nbs.vec[p*nnbs+i] = neighbors.vec[cum_n_nbs+n_nbs-1];
        }
    }

    nbs.copy_vec_to_ptr ();

}


void Points :: write (const Io& io) const
{
    Double2 position_buffer (npoints, Double1(3));
    Double2 velocity_buffer (npoints, Double1(3));
//
    for (size_t p = 0; p < npoints; p++)
    {
        position_buffer[p] = {position.vec[p].x(),
                              position.vec[p].y(),
                              position.vec[p].z() };
        velocity_buffer[p] = {velocity.vec[p].x(),
                              velocity.vec[p].y(),
                              velocity.vec[p].z() };
    }

    io.write_array (prefix+"position", position_buffer);
    io.write_array (prefix+"velocity", velocity_buffer);

//    Long1 n_neighbors_buffer (npoints);
//    Long1   neighbors_buffer (tot_n_neighbors);

//    for (size_t p = 0; p < npoints; p++)
//    {
//        n_neighbors_buffer[p] = n_neighbors[p];
//    }

//    for (size_t p = 0; p < tot_n_neighbors; p++)
//    {
//        neighbors_buffer[p] = neighbors[p];
//    }

//    io.write_list (prefix+"n_neighbors", n_neighbors_buffer);
//    io.write_list (prefix+  "neighbors",   neighbors_buffer);

    io.write_list (prefix+"n_neighbors", n_neighbors.vec);
    io.write_list (prefix+  "neighbors",   neighbors.vec);
}
