#include "model.hpp"


void Model :: read (const Io& io)
{
//    chemistry.read (io);
     geometry.read (io);

//     lines.set_npoints (geometry.get_npoints());

//        lines.read (io);

    thermodynamics.read (io);

}


void Model :: write (const Io& io) const
{
//    chemistry.write (io);
     geometry.write (io);
//        lines.write (io);

    thermodynamics.write (io);
}