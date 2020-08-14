#include "rays.hpp"


const string prefix = "geometry/rays/";


void Rays :: read (const Io& io)
{
    cout << "Reading rays..." << endl;

    io.read_length (prefix+"direction", nrays);

    direction = (Vector*) pc::accelerator::malloc (nrays*sizeof(Vector));
    antipod   = (size_t*) pc::accelerator::malloc (nrays*sizeof(size_t));
    weight    = (double*) pc::accelerator::malloc (nrays*sizeof(double));

    Double2 direction_buffer (nrays, Double1(3));
    Double1    weight_buffer (nrays);

    io.read_array (prefix+"direction", direction_buffer);
    io.read_list  (prefix+"weight",       weight_buffer);

    for (size_t r = 0; r < nrays; r++)
    {
        direction[r] = Vector (direction_buffer[r][0],
                               direction_buffer[r][1],
                               direction_buffer[r][2] );

           weight[r] = weight_buffer[r];
    }

    cout << "nrays = " << nrays << endl;

    const double tolerance = 1.0E-9;

    for (size_t r1 = 0; r1 < nrays; r1++)
    {
        for (size_t r2 = 0; r2 < nrays; r2++)
        {
            if ((direction[r1] + direction[r2]).squaredNorm() < tolerance)
            {
                antipod[r1] = r2;
            }
        }
    }
}


void Rays :: write (const Io& io) const
{
    Double2 direction_buffer (nrays, Double1(3));
    Double1    weight_buffer (nrays);

    for (size_t r = r; r < nrays; r++)
    {
        direction_buffer[r] = {direction[r].x(),
                               direction[r].y(),
                               direction[r].z() };

           weight_buffer[r] = weight[r];
    }

    io.read_array (prefix+"direction", direction_buffer);
    io.read_list  (prefix+"weight",       weight_buffer);
}
