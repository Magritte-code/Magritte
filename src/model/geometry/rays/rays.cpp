#include "rays.hpp"


const string prefix = "geometry/rays/";


void Rays :: read (const Io& io)
{
    cout << "Reading rays..." << endl;

    parameters->set_nrays (io.get_length (prefix+"direction"));
    parameters->set_hnrays (parameters->nrays()/2);

    direction.resize (parameters->nrays());
      antipod.resize (parameters->nrays());
       weight.resize (parameters->nrays());

    Real2 direction_buffer (parameters->nrays(), Real1(3));

    io.read_array (prefix+"direction", direction_buffer);
    io.read_list  (prefix+"weight",    weight);

    for (Size r = 0; r < parameters->nrays(); r++)
    {
        direction[r] = Vector3D (direction_buffer[r][0],
                                 direction_buffer[r][1],
                                 direction_buffer[r][2] );
    }

    const Real tolerance = 1.0E-9;

    for (Size r1 = 0; r1 < parameters->nrays(); r1++)
    {
        for (Size r2 = 0; r2 < parameters->nrays(); r2++)
        {
            if ((direction[r1] + direction[r2]).squaredNorm() < tolerance)
            {
                antipod[r1] = r2;
            }
        }
    }

    direction.copy_vec_to_ptr ();
    antipod  .copy_vec_to_ptr ();
    weight   .copy_vec_to_ptr ();
}


void Rays :: write (const Io& io) const
{
    cout << "Writing rays..." << endl;

    Real2 direction_buffer (parameters->nrays(), Real1(3));

    for (Size r = 0; r < parameters->nrays(); r++)
    {
        direction_buffer[r] = {direction[r].x(),
                               direction[r].y(),
                               direction[r].z() };
    }

    io.write_array (prefix+"direction", direction_buffer);
    io.write_list  (prefix+"weight",    weight          );
}
