#include "geometry.hpp"


void Geometry :: read (const Io& io)
{
    points.read (io);

    npoints = points.get_npoints();
    boundary.set_npoints (npoints);

    cout << npoints << endl;

    rays.read (io);

    nrays = rays.get_nrays();

    boundary.read (io);

    cout << "Resizing lengths" << endl;

    lengths.resize (nrays/2*npoints);

    cout << "Done Reading!" << endl;
}


void Geometry :: write (const Io& io) const
{
    points  .write (io);
    rays    .write (io);
    boundary.write (io);
}
