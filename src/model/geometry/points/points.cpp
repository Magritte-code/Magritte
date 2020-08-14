#include "points.hpp"


const string prefix = "geometry/points/";


void Points :: read (const Io& io)
{
    cout << "Reading points..." << endl;

    size_t n1, n2;
    io.read_length (prefix+"position", n1);
    io.read_length (prefix+"velocity", n2);
    assert (n1 == n2);
    npoints = n1;

    cout << "npoints = " << npoints << endl;

    position = (Vector*) pc::accelerator::malloc (npoints*sizeof(Vector));
    velocity = (Vector*) pc::accelerator::malloc (npoints*sizeof(Vector));

    Double2 position_buffer (npoints, Double1(3));
    Double2 velocity_buffer (npoints, Double1(3));

    io.read_array (prefix+"position", position_buffer);
    io.read_array (prefix+"velocity", velocity_buffer);

    for (size_t p = 0; p < npoints; p++)
    {
        position[p] = Vector (position_buffer[p][0],
                              position_buffer[p][1],
                              position_buffer[p][2] );
        velocity[p] = Vector (velocity_buffer[p][0],
                              velocity_buffer[p][1],
                              velocity_buffer[p][2] );
    }

    io.read_length (prefix+"neighbors", tot_n_neighbors);

    cout << "tot_n_neighbors = " << tot_n_neighbors << endl;

    cum_n_neighbors = (size_t*) pc::accelerator::malloc (npoints        *sizeof(size_t));
        n_neighbors = (size_t*) pc::accelerator::malloc (npoints        *sizeof(size_t));
          neighbors = (size_t*) pc::accelerator::malloc (tot_n_neighbors*sizeof(size_t));

    cout << "memory_allocated = " << endl;

    Long1 n_neighbors_buffer (npoints);
    Long1   neighbors_buffer (tot_n_neighbors);

    cout << "lists made = " << endl;

    io.read_list (prefix+"n_neighbors", n_neighbors_buffer);
    io.read_list (prefix+  "neighbors",   neighbors_buffer);

    cout << "lists read = " << endl;

        n_neighbors[0] = n_neighbors_buffer[0];
    cum_n_neighbors[0] = 0;

    cout << "first put" << endl;

    for (size_t p = 1; p < npoints; p++)
    {
            n_neighbors[p] = n_neighbors_buffer[p];
        cum_n_neighbors[p] = cum_n_neighbors[p-1] + n_neighbors[p-1];
    }

    cout << "points put" << endl;
    cout << "n " << neighbors_buffer.size() << endl;

    for (size_t p = 0; p < tot_n_neighbors; p++)
    {
//        cout << p << endl;
        neighbors[p] = neighbors_buffer[p];
    }

    cout << "neighbors put" << endl;
}


void Points :: write (const Io& io) const
{
    Double2 position_buffer (npoints, Double1(3));
    Double2 velocity_buffer (npoints, Double1(3));

    for (size_t p = 0; p < npoints; p++)
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

    Long1 n_neighbors_buffer (npoints);
    Long1   neighbors_buffer (tot_n_neighbors);

    for (size_t p = 0; p < npoints; p++)
    {
        n_neighbors_buffer[p] = n_neighbors[p];
    }

    for (size_t p = 0; p < tot_n_neighbors; p++)
    {
        neighbors_buffer[p] = neighbors[p];
    }

    io.write_list (prefix+"n_neighbors", n_neighbors_buffer);
    io.write_list (prefix+  "neighbors",   neighbors_buffer);
}
