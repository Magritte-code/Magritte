#include "image.hpp"
#include "paracabs.hpp"


const string prefix = "image/";


///  Constructor for Image
//////////////////////////
Image :: Image (const Geometry& geometry, const Size rr) : ray_nr (rr)
{
    if (geometry.parameters.dimension() == 1)
    {
        if ((geometry.rays.direction[ray_nr].x() != 0.0) ||
            (geometry.rays.direction[ray_nr].y() != 1.0) ||
            (geometry.rays.direction[ray_nr].z() != 0.0)   )
        {
            throw std::runtime_error ("In 1D, the image ray has to be (0,1,0)");
        }
    }

    ImX.resize (geometry.parameters.npoints());
    ImY.resize (geometry.parameters.npoints());
    I.  resize (geometry.parameters.npoints(), geometry.parameters.nfreqs());

    set_coordinates (geometry);
}


///  Copy constructor for Image
///////////////////////////////
Image :: Image (const Image& image) : ray_nr (image.ray_nr)
{
    ImX = image.ImX;
    ImY = image.ImY;

    // Deep copy of I
    I.nrows = image.I.nrows;
    I.ncols = image.I.ncols;

    I.vec            = image.I.vec;
    I.allocated      = false;
    I.allocated_size = 0;
    I.set_dat ();
}

///  print: write out the images
///    @param[in] io: io object
////////////////////////////////
//void Image :: write (const Io &io) const
//{
//    cout << "Writing image..."    << endl;
//
//    const string str_ray_nr = std::to_string (ray_nr);
//
//    io.write_list  (prefix+"ImX_"+str_ray_nr, ImX);
//    io.write_list  (prefix+"ImY_"+str_ray_nr, ImY);
//
//    // Create one intensity variable and reuse for I_m and I_p to save memory.
//    Double2 intensity_m (ncells, Double1 (nfreqs));
//    Double2 intensity_p (ncells, Double1 (nfreqs));
//
//    OMP_PARALLEL_FOR (p, ncells)
//    {
//        for (size_t f = 0; f < nfreqs; f++)
//        {
//            intensity_m[p][f] = get (I_m[p], f);
//            intensity_p[p][f] = get (I_p[p], f);
//        }
//    }
//
//    io.write_array (prefix+"I_m_"+str_ray_nr, intensity_m);
//    io.write_array (prefix+"I_p_"+str_ray_nr, intensity_p);
//}


///  Setter for the coordinates on the image axes
///    @param[in] geometry : geometry object of the model
/////////////////////////////////////////////////////////
void Image :: set_coordinates (const Geometry& geometry)
{
    if (geometry.parameters.dimension() == 1)
    {
        threaded_for (p, geometry.parameters.npoints(),
        {
            ImX[p] = geometry.points.position[p].x();
            ImY[p] = 0.0;
        })
    }

    if (geometry.parameters.dimension() == 3)
    {
        const double rx = geometry.rays.direction[ray_nr].x();
        const double ry = geometry.rays.direction[ray_nr].y();
        const double rz = geometry.rays.direction[ray_nr].z();

        const double         denominator = sqrt (rx*rx + ry*ry);
        const double inverse_denominator = 1.0 / denominator;

        const double ix =  ry * inverse_denominator;
        const double iy = -rx * inverse_denominator;

        const double jx =  rx * rz * inverse_denominator;
        const double jy =  ry * rz * inverse_denominator;
        const double jz = -denominator;

        if (denominator != 0.0)
        {
            threaded_for (p, geometry.parameters.npoints(),
            {
                ImX[p] =   ix * geometry.points.position[p].x()
                         + iy * geometry.points.position[p].y();

                ImY[p] =   jx * geometry.points.position[p].x()
                         + jy * geometry.points.position[p].y()
                         + jz * geometry.points.position[p].z();

                // cout << "position = " << geometry.points.position[p].x() << endl;
                // cout << "ImX = " << ImX[p] << "  " << ix << endl;
                // cout << "ImY = " << ImY[p] << "  " << jx << endl;
            })
        }
        else
        {
            threaded_for (p, geometry.parameters.npoints(),
            {
                ImX[p] = geometry.points.position[p].x();
                ImY[p] = geometry.points.position[p].y();
            })
        }
    }
}
