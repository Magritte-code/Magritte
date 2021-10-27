#include "rays.hpp"


const string prefix = "geometry/rays/";


void Rays :: read (const Io& io)
{
    cout << "Reading rays..." << endl;

    parameters.set_nrays (io.get_length (prefix+"direction"));
    parameters.set_hnrays (parameters.nrays()/2);

    direction.resize (parameters.nrays());
    weight   .resize (parameters.nrays());
    antipod  .resize (parameters.nrays());

    Double2 direction_buffer (parameters.nrays(), Double1(3));

    io.read_array (prefix+"direction", direction_buffer);
    io.read_list  (prefix+"weight",    weight);

    for (Size r = 0; r < parameters.nrays(); r++)
    {
        direction[r] = Vector3D (direction_buffer[r][0],
                                 direction_buffer[r][1],
                                 direction_buffer[r][2] );
    }

    const double tolerance = 1.0E-9;

    for (Size r1 = 0; r1 < parameters.nrays(); r1++)
    {
        for (Size r2 = 0; r2 < parameters.nrays(); r2++)
        {
            if ((direction[r1] + direction[r2]).squaredNorm() < tolerance)
            {
                antipod[r1] = r2;
            }
        }
    }

    direction.copy_vec_to_ptr ();
    weight   .copy_vec_to_ptr ();
    antipod  .copy_vec_to_ptr ();
}


void Rays :: write (const Io& io) const
{
    cout << "Writing rays..." << endl;

    Double2 direction_buffer (parameters.nrays(), Double1(3));

    for (Size r = 0; r < parameters.nrays(); r++)
    {
        direction_buffer[r] = {direction[r].x(),
                               direction[r].y(),
                               direction[r].z() };
    }

    io.write_array (prefix+"direction", direction_buffer);
    io.write_list  (prefix+"weight",    weight          );
}


void Rays :: set_1D_adaptive_rays (const Vector<Vector3D>& position)
{
    // Convenience parameter
    const Size N = parameters.npoints();

    // Set nrays and half nrays
    parameters.set_nrays  (2*(N+2));
    parameters.set_hnrays (parameters.nrays()/2);

    // Allocate memory
    m_direction.resize (N, parameters.nrays());
    m_weight   .resize (N, parameters.nrays());
    antipod    .resize (   parameters.nrays());


    // For all different radii
    threaded_for (o, N,
    {
        const Real p_o = position[o].x();


        // Set ray directions
        /////////////////////

        // First ray
        m_direction(o,0) = Vector3D (1.0, 0.0, 0.0);

        // Middle rays
        for (Size r = 0; r < N; r++)
        {
            const Real p_r = position[r].x();

            const Real rx = p_o / sqrt(p_o*p_o + p_r*p_r);
            const Real ry = p_r / sqrt(p_o*p_o + p_r*p_r);

            m_direction(o,r+1) = Vector3D (rx, ry, 0.0);
        }

        // Last ray
        m_direction(o,N+1) = Vector3D (0.0, 1.0, 0.0);


        // Set ray weights
        //////////////////
    
        // First ray
        Vector3D upper = m_direction(o,0) + m_direction(o,1);
        Vector3D lower = m_direction(o,0);

        m_weight(o,0) =   0.5*upper.y() / sqrt(upper.squaredNorm())\
                        - 0.5*lower.y() / sqrt(lower.squaredNorm());


        // Middle rays
        for (Size r = 1; r < N+1; r++)
        {
            upper = m_direction(o,r) + m_direction(o,r+1);
            lower = m_direction(o,r) + m_direction(o,r-1);

            m_weight(o,r) =   0.5*upper.y() / sqrt(upper.squaredNorm())\
                            - 0.5*lower.y() / sqrt(lower.squaredNorm());
        }

        // Last ray
        upper = m_direction(o,N+1);
        lower = m_direction(o,N+1) + m_direction(o,N);

        m_weight(o,N+1) =   0.5*upper.y() / sqrt(upper.squaredNorm())\
                          - 0.5*lower.y() / sqrt(lower.squaredNorm());
        

        // Set antipodal rays
        /////////////////////
        
        for (Size r = 0; r < N+2; r++)
        {
            const Real rx = m_direction(o,r).x();
            const Real ry = m_direction(o,r).y();

            m_direction(o,N+2+r) = Vector3D(-rx,-ry,0.0);
            m_weight   (o,N+2+r) = m_weight(o,r);

            antipod[N+2+r] = r;
            antipod[r] = N+2+r;
        }
    })


    m_direction.copy_vec_to_ptr ();
    m_weight   .copy_vec_to_ptr ();
    antipod    .copy_vec_to_ptr ();
}
