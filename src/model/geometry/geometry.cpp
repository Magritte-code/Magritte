#include "geometry.hpp"


void Geometry :: read (const Io& io)
{
    points  .read (io);
    rays    .read (io);
    boundary.read (io);

    cout << "Resizing lengths" << endl;
    cout << parameters.hnrays() * parameters.npoints() << endl;

    lengths.resize (parameters.hnrays()*parameters.npoints());

    cout << "Done Reading!" << endl;
}


void Geometry :: write (const Io& io) const
{
    points  .write (io);
    rays    .write (io);
    boundary.write (io);
}
