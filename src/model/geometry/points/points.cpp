#include <assert.h>
#include "points.hpp"
#include "tools/types.hpp"


const string prefix = "geometry/points/";


void Points :: read (const Io& io)
{
    cout << "Reading points..." << endl;

    parameters.set_npoints (io.get_length (prefix+"position"));
    parameters.set_npoints (io.get_length (prefix+"velocity"));

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
    vector <Size> n_neighbors;
    vector <Size>   neighbors;

    n_neighbors.resize (parameters.npoints());
      neighbors.resize (parameters.totnnbs());

    cout << "memory_allocated = " << endl;

    cout << "lists made = " << endl;

    cout << &n_neighbors << endl;


    io.read_list (prefix+"n_neighbors", n_neighbors);
    io.read_list (prefix+  "neighbors",   neighbors);

    cout << n_neighbors[0] << endl;
    cout << neighbors[0] << endl;

    this->multiscale.set_all_neighbors(n_neighbors,neighbors);


    cout << "lists read = " << endl;


    cout << "first put" << endl;

    cout << &n_neighbors << endl;
    cout << &n_neighbors[0]        << endl;
    cout << &n_neighbors      << endl;


    position.copy_vec_to_ptr ();
    velocity.copy_vec_to_ptr ();

}


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

    io.write_list (prefix+"n_neighbors", this->multiscale.get_all_n_neighbors());
    io.write_list (prefix+  "neighbors", this->multiscale.get_all_neighbors_as_vector());

}
