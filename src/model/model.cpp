#include "model.hpp"


void Model :: read (const Io& io)
{
    geometry.read (io);
}


void Model :: write (const Io& io) const
{
    geometry.write (io);
}