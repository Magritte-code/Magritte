#include <assert.h>
#include "points.hpp"


const string prefix = "geometry/points/";


void Points :: read (const Io& io)
{
    cout << "Reading points..." << endl;

/*@Frederik: either these two try to set it to the same value, in which case nothing happens;
or they try to set it to a different value, in which case a DoubleSetException is thrown, so i'm commenting one out for now
Also assert that the lengths are the same; You can change this to any other exception
*/
    assert(io.get_length (prefix+"position")==io.get_length (prefix+"velocity"));
    parameters.set_npoints (io.get_length (prefix+"position"));
//    parameters.set_npoints (io.get_length (prefix+"velocity"));

    cout << "npoints = " << parameters.npoints() << endl;

    position.resize (parameters.npoints());
    velocity.resize (parameters.npoints());

    Double2 position_buffer (parameters.npoints(), Double1(3));
    Double2 velocity_buffer (parameters.npoints(), Double1(3));

    io.read_array (prefix+"position", position_buffer);
    io.read_array (prefix+"velocity", velocity_buffer);

    for (Size p = 0; p < parameters.npoints(); p++)
    {
        position[p] = Vector3D (position_buffer[p][0],
                                position_buffer[p][1],
                                position_buffer[p][2] );

        velocity[p] = Vector3D (velocity_buffer[p][0],
                                velocity_buffer[p][1],
                                velocity_buffer[p][2] );
    }

    parameters.set_totnnbs (io.get_length (prefix+"neighbors"));

    cout << "tot_n_neighbors = " << parameters.totnnbs() << endl;

//Temporary vectors for neighbors
    Vector <Size> n_neighbors;
    Vector <Size>   neighbors;

//    cum_n_neighbors.resize (parameters.npoints());
    n_neighbors.resize (parameters.npoints());
      neighbors.resize (parameters.totnnbs());

    cout << "memory_allocated = " << endl;

    cout << "lists made = " << endl;

    cout << n_neighbors.vec.data() << endl;

    io.read_list (prefix+"n_neighbors", n_neighbors);
    io.read_list (prefix+  "neighbors",   neighbors);

    self.curr_neighbors=Neighbors(n_neighbors,neighbors);

    cout << "lists read = " << endl;

//    cum_n_neighbors[0] = 0;

    cout << "first put" << endl;

    cout << n_neighbors.vec.data() << endl;
    cout << &n_neighbors[0]        << endl;
    cout << n_neighbors.dat       << endl;
    cout << n_neighbors.ptr       << endl;

//    for (Size p = 1; p < parameters.npoints(); p++)
//    {
//        cum_n_neighbors[p] = cum_n_neighbors[p-1] + n_neighbors[p-1];
//    }

    cout << "points put" << endl;

    position.copy_vec_to_ptr ();
    velocity.copy_vec_to_ptr ();

//    cum_n_neighbors.copy_vec_to_ptr ();
//        n_neighbors.copy_vec_to_ptr ();
//          neighbors.copy_vec_to_ptr ();

    cout << "neighbors put" << endl;

//@Frederik: nbs does not seem to do anything remotely useful
//     nbs.resize (parameters.npoints()*nnbs);
//
//     for (Size p = 0; p < parameters.npoints(); p++)
//     {
//         const Size     n_nbs =     n_neighbors[p];
//         const Size cum_n_nbs = cum_n_neighbors[p];
//
//         for (Size i = 0; (i < n_nbs) && (i < nnbs); i++)
//         {
//             nbs[p*nnbs+i] = neighbors[cum_n_nbs+i];
//         }
//         for (Size i = n_nbs; i < nnbs; i++)
//         {
//             nbs[p*nnbs+i] = neighbors[cum_n_nbs+n_nbs-1];
//         }
//     }
//
//     nbs.copy_vec_to_ptr ();
// }


void Points :: write (const Io& io) const
{
    Double2 position_buffer (parameters.npoints(), Double1(3));
    Double2 velocity_buffer (parameters.npoints(), Double1(3));

    for (size_t p = 0; p < parameters.npoints(); p++)
    {
        position_buffer[p] = {position[p].x(),
                              position[p].y(),
                              position[p].z() };
        velocity_buffer[p] = {velocity[p].x(),
                              velocity[p].y(),
                              velocity[p].z() };
    }

    io.write_array (prefix+"position", position_buffer);
    io.write_array (prefix+"velocity", velocity_buffer);

//@Frederik: this might write different results depending on the current level of coarsening
//TODO take another look at this
    io.write_list (prefix+"n_neighbors", self.curr_neighbors.n_neighbors);
    io.write_list (prefix+  "neighbors", self.curr_neighbors.get_flattened_neigbors_list());
}
