#include "geometry.hpp"


void Geometry :: read (const Io& io)
{
    points.read (io);

    if (parameters.adaptive_ray_tracing())
    {
        if ((parameters.dimension() == 1) && (parameters.spherical_symmetry()))
        {
            rays.set_1D_adaptive_rays (points.position);
        }
        else
        {
            throw std::runtime_error ("Adaptive ray-tracing is only implemented in 1D.");
        }
    }
    else
    {
        rays.read (io);
    }

    boundary.read (io);

    lengths.resize (parameters.hnrays(), parameters.npoints());
}


void Geometry :: write (const Io& io) const
{
    points.write (io);

    if (!parameters.adaptive_ray_tracing())
    {
        rays.write (io);
    }

    boundary.write (io);
}
