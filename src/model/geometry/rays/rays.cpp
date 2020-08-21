#include "rays.hpp"


const string prefix = "geometry/rays/";


void Rays :: read (const Io& io)
{
    cout << "Reading rays..." << endl;

    io.read_length (prefix+"direction", nrays);

    direction.resize (nrays);
      antipod.resize (nrays);
       weight.resize (nrays);

//    direction = (Vector3D*) pc::accelerator::malloc (nrays*sizeof(Vector3D));
//    antipod   = (size_t*)   pc::accelerator::malloc (nrays*sizeof(size_t));
//    weight    = (double*)   pc::accelerator::malloc (nrays*sizeof(double));

    Double2 direction_buffer (nrays, Double1(3));
//    Double1    weight_buffer (nrays);

    io.read_array (prefix+"direction", direction_buffer);
//    io.read_list  (prefix+"weight",       weight_buffer);
    io.read_list  (prefix+"weight",    weight.vec);

    for (size_t r = 0; r < nrays; r++)
    {
        direction.vec[r] = Vector3D (direction_buffer[r][0],
                                     direction_buffer[r][1],
                                     direction_buffer[r][2] );
//
//           weight[r] = weight_buffer[r];
    }

    cout << "nrays = " << nrays << endl;

    const double tolerance = 1.0E-9;

    for (size_t r1 = 0; r1 < nrays; r1++)
    {
        for (size_t r2 = 0; r2 < nrays; r2++)
        {
            if ((direction.vec[r1] + direction.vec[r2]).squaredNorm() < tolerance)
            {
                antipod.vec[r1] = r2;
            }
        }
    }

    direction.copy_vec_to_ptr ();
    antipod  .copy_vec_to_ptr ();
    weight   .copy_vec_to_ptr ();
}


void Rays :: write (const Io& io) const
{
    Double2 direction_buffer (nrays, Double1(3));
//    Double1    weight_buffer (nrays);

    for (Size r = 0; r < nrays; r++)
    {
        direction_buffer[r] = {direction.vec[r].x(),
                               direction.vec[r].y(),
                               direction.vec[r].z() };

//           weight_buffer[r] = weight[r];
    }

    io.write_array (prefix+"direction", direction_buffer);
//    io.write_list  (prefix+"weight",       weight_buffer);
    io.write_list  (prefix+"weight",    weight.vec);
}
