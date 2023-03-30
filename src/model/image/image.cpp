#include "image.hpp"
#include "paracabs.hpp"
#include <limits>
#include <algorithm>


const string prefix = "image/";


///  Constructor for Image; using default rays, default point position
//////////////////////////
Image :: Image (const Geometry& geometry, const ImageType it, const Size rr) : imageType(it), imagePointPosition(AllModelPoints), ray_nr(rr), ray_direction(Vector3D(geometry.rays.direction[ray_nr].x(), geometry.rays.direction[ray_nr].y(), geometry.rays.direction[ray_nr].z()))
{
    if (geometry.parameters->dimension() == 1)
    {
        if ((geometry.rays.direction[ray_nr].x() != 0.0) ||
            (geometry.rays.direction[ray_nr].y() != 1.0) ||
            (geometry.rays.direction[ray_nr].z() != 0.0)   )
        {
            throw std::runtime_error ("In 1D, the image ray has to be (0,1,0)");
        }
    }

    // ImX.resize (geometry.parameters->npoints());
    // ImY.resize (geometry.parameters->npoints());
    // I.  resize (geometry.parameters->npoints(), geometry.parameters->nfreqs());
    // ray_origin.resize(geometry.parameters->npoints());
    // ray_direction = Vector3D(geometry.rays.direction[ray_nr].x(), geometry.rays.direction[ray_nr].y(), geometry.rays.direction[ray_nr].z());

    set_coordinates_all_model_points (geometry);
}

///  Constructor for Image; using default rays, non-default point position
//////////////////////////
Image :: Image (const Geometry& geometry, const ImageType it, const Size rr, const Size Nxpix, const Size Nypix) : imageType(it), imagePointPosition(ProjectionSurface), ray_nr(rr), ray_direction(Vector3D(geometry.rays.direction[ray_nr].x(), geometry.rays.direction[ray_nr].y(), geometry.rays.direction[ray_nr].z()))
{
    if (geometry.parameters->dimension() == 1)
    {
        if ((geometry.rays.direction[ray_nr].x() != 0.0) ||
            (geometry.rays.direction[ray_nr].y() != 1.0) ||
            (geometry.rays.direction[ray_nr].z() != 0.0)   )
        {
            throw std::runtime_error ("In 1D, the image ray has to be (0,1,0)");
        }
    }

    // ImX.resize (geometry.parameters->npoints());
    // ImY.resize (geometry.parameters->npoints());
    // I.  resize (geometry.parameters->npoints(), geometry.parameters->nfreqs());
    // ray_origin.resize(geometry.parameters->npoints());
    // ray_direction = const Vector3D(geometry.rays.direction[ray_nr].x(), geometry.rays.direction[ray_nr].y(), geometry.rays.direction[ray_nr].z());

    set_coordinates_projection_surface (geometry, Nxpix, Nypix);
}


///  Copy constructor for Image
///////////////////////////////
Image :: Image (const Image& image) : imageType(image.imageType), imagePointPosition(image.imagePointPosition), ray_nr (image.ray_nr), ray_direction(image.ray_direction)
{
    ImX = image.ImX;
    ImY = image.ImY;
    //FIX ME: check whether this copy constructor actually properly copies these values!
    // ray_origin = image.ray_origin;
    // ray_direction = image.ray_direction;


    // Deep copy of I
    I.nrows = image.I.nrows;
    I.ncols = image.I.ncols;

    I.vec            = image.I.vec;
    I.allocated      = false;
    I.allocated_size = 0;
    I.set_dat ();
}

///  Constructor for Image; using non-default rays, default point position
//////////////////////////
Image :: Image (const Geometry& geometry, const ImageType it, const Vector3D raydir) : imageType(it), imagePointPosition(AllModelPoints), ray_nr(-1), ray_direction(raydir), closest_bdy_point(geometry.parameters->npoints())
{
    if (geometry.parameters->dimension() == 1)
    {
        // throw std::runtime_error ("1D is currently not supported, when choosing a custom ray direction.");
        //TODO: either do a lot of work (do everything correctly in )
        //OR: set the raydirection just to [0,1,0]
        //Dev note: in Geometry, I should define how to handle custom ray directions in 1D...
        if ((raydir.x() != 0.0) ||
            (raydir.y() != 1.0) ||
            (raydir.z() != 0.0)   )
        {
            throw std::runtime_error ("In 1D, the image ray has to be (0,1,0)");
        }
    }

    // ImX.resize (geometry.parameters->npoints());
    // ImY.resize (geometry.parameters->npoints());
    // I.  resize (geometry.parameters->npoints(), geometry.parameters->nfreqs());
    // ray_origin.resize(geometry.parameters->npoints());

    // ray_direction = raydir;

    set_coordinates_all_model_points (geometry);
}


///  Constructor for Image; using non-default rays, non-default point position
//////////////////////////
Image :: Image (const Geometry& geometry, const ImageType it, const Vector3D raydir, const Size Nxpix, const Size Nypix) : imageType(it), imagePointPosition(ProjectionSurface), ray_nr(-1), ray_direction(raydir)
{
    if (geometry.parameters->dimension() == 1)
    {
        // throw std::runtime_error ("1D is currently not supported, when choosing a custom ray direction.");
        //TODO: either do a lot of work (do everything correctly in )
        //OR: set the raydirection just to [0,1,0]
        //Dev note: in Geometry, I should define how to handle custom ray directions in 1D...
        if ((raydir.x() != 0.0) ||
            (raydir.y() != 1.0) ||
            (raydir.z() != 0.0)   )
        {
            throw std::runtime_error ("In 1D, the image ray has to be (0,1,0)");
        }
    }

    // ImX.resize (geometry.parameters->npoints());
    // ImY.resize (geometry.parameters->npoints());
    // I.  resize (geometry.parameters->npoints(), geometry.parameters->nfreqs());
    // ray_origin.resize(geometry.parameters->npoints());

    // ray_direction = raydir;

    set_coordinates_projection_surface (geometry, Nxpix, Nypix);
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
inline void Image :: set_coordinates_all_model_points (const Geometry& geometry)
{
    ImX.resize (geometry.parameters->npoints());
    ImY.resize (geometry.parameters->npoints());
    I.  resize (geometry.parameters->npoints(), geometry.parameters->nfreqs());
    // ray_origin.resize(geometry.parameters->npoints());
    closest_bdy_point = geometry.parameters->npoints();

    if (geometry.parameters->dimension() == 1)
    {
        threaded_for (p, geometry.parameters->npoints(),
        {
            ImX[p] = geometry.points.position[p].x();
            ImY[p] = 0.0;
        })
    }

    if (geometry.parameters->dimension() == 3)
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

        if (denominator >= 1.0e-9)
        {
            threaded_for (p, geometry.parameters->npoints(),
            {
                ImX[p] =   ix * geometry.points.position[p].x()
                         + iy * geometry.points.position[p].y();

                ImY[p] =   jx * geometry.points.position[p].x()
                         + jy * geometry.points.position[p].y()
                         + jz * geometry.points.position[p].z();
            })
        }
        else
        {
            threaded_for (p, geometry.parameters->npoints(),
            {
                ImX[p] = geometry.points.position[p].x();
                ImY[p] = geometry.points.position[p].y();
            })
        }
    }
}

inline void Image :: set_coordinates_projection_surface (const Geometry& geometry, const Size Nxpix, const Size Nypix)//(const Geometry& geometry)
{
    //If the raydirection is [0,0,1], then I set hat(x)=[1,0,0], hat(y)=[0,1,0] (lefthand rule: x, n, y)
    //For other ray directions, I define the (almost everywhere) unique rotation as follows:
    //First, we rotate around the (default) y vector (with angle α), then around the x vector (with angle β)
    //This results in the following rotation matrix:
    //[cos(α)       , 0     , -sin(α)      ]
    //[-sin(α)sin(β), cos(β), -cos(α)sin(β)]
    //[sin(α)cos(β) , sin(β), cos(α)cos(β) ]
    //Given the raydirection (xn,yn,zn), the rotation angles can be found by looking at the last column
    //α = arcsin(-xn)
    //if x!=+-1 (i.e. √(1-xn^2)>0)
    //yn = -cos(arcsin(-xn))*sin(β) = -√(1-xn^2)*sin(β)
    //so sin(β) = -yn/√(1-xn^2)
    //or cos(β) = zn/√(1-xn^2)
    //
    //Translating this to the hat(x) and hat(y) vectors:
    //for hat(x):
    //xx = cos(α) = cos(arcsin(-xn)) = √(1-xn^2)
    //yx = -sin(α)sin(β) = xn * -yn/√(1-xn^2)
    //zx = sin(α)cos(β) = -xn * zn/√(1-xn^2)
    //for hat(y):
    //xy = 0
    //yy = cos(β) = zn/√(1-xn^2)
    //yz = sin(β) = -yn/√(1-xn^2)

    //OR WE USE THE DIRECTIONS DEFINED ABOVE (only difference is how the axes were defined):

    //TODO copy


    //Then figure out (using the boundary points, which points lie farthest in x-y space, then use to define rectangular grid? (at least min/max coords))

    //default initialization for the limits of the projected boundary points
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    //Maybe in other function: define the imaging plane (do we still need the perpendicular vectors anywhere? I guess not)
    if (geometry.parameters->dimension() == 1)
    {//TODO: figure out what we actually want in 1D!
        for (Size bdy_idx = 0; bdy_idx<geometry.parameters->nboundary(); bdy_idx++)
        {
            const Size bdy_point_index = geometry.boundary.boundary2point[bdy_idx];
            const double ImX = geometry.points.position[bdy_point_index].x();
            max_x = std::max(max_x, ImX);
            min_x = std::min(min_x, ImX);
        }
        // threaded_for (p, geometry.parameters->npoints(),
        // {
        //     ImX[p] = geometry.points.position[p].x();
        //     // ImY[p] = 0.0;
        // })
        max_x = std::max(std::abs(min_x), max_x);
        min_x = -max_x;
        min_y = min_x;
        max_y = max_x;
    }

    if (geometry.parameters->dimension() == 3)
    {
        const double rx = ray_direction.x();
        const double ry = ray_direction.y();
        const double rz = ray_direction.z();
        std::cout<<"raydir: "<<rx<<" "<<ry<<" "<<rz<<std::endl;

        const double         denominator = sqrt (rx*rx + ry*ry);
        const double inverse_denominator = 1.0 / denominator;

        const double ix =  ry * inverse_denominator;
        const double iy = -rx * inverse_denominator;

        const double jx =  rx * rz * inverse_denominator;
        const double jy =  ry * rz * inverse_denominator;
        const double jz = -denominator;



        if (denominator >= 1.0e-9)
        {
            //Err, just find the projected bounds of the boundary points
            //TODO: try to find paracabs fun that actually implements reduce operation (or use the ThreadPrivate stuff, then manually min, max them)
            for (Size bdy_idx = 0; bdy_idx<geometry.parameters->nboundary(); bdy_idx++)
            {
                const Size bdy_point_index = geometry.boundary.boundary2point[bdy_idx];
                const double ImX = ix * geometry.points.position[bdy_point_index].x()
                                 + iy * geometry.points.position[bdy_point_index].y();
                const double ImY = jx * geometry.points.position[bdy_point_index].x()
                                 + jy * geometry.points.position[bdy_point_index].y()
                                 + jz * geometry.points.position[bdy_point_index].z();
                min_x = std::min(min_x, ImX);
                max_x = std::max(max_x, ImX);
                min_y = std::min(min_y, ImY);
                max_y = std::max(max_y, ImY);
            }
        }
        else
        {
            for (Size bdy_idx = 0; bdy_idx<geometry.parameters->nboundary(); bdy_idx++)
            {
                const Size bdy_point_index = geometry.boundary.boundary2point[bdy_idx];
                const double ImX = geometry.points.position[bdy_point_index].x();
                const double ImY = geometry.points.position[bdy_point_index].y();

                max_x = std::max(max_x, ImX);
                max_y = std::max(max_y, ImY);
                min_x = std::min(min_x, ImX);
                min_y = std::max(min_y, ImY);
            }
        }
    }

    //Define pixels spanned by the plane: TODO: figure out why I thought I needed the closest bdy point
    closest_bdy_point = geometry.get_closest_bdy_point_in_custom_raydir(ray_direction);
    const Vector3D closest_bdy_position = geometry.points.position[closest_bdy_point];
    double distance_from_origin = closest_bdy_position.dot(ray_direction);//in a specific ray direction
    //for 1D spherical symmetry, this is incorrect, as we represent the distance from the middle by the x-coordinate
    if (geometry.parameters->spherical_symmetry())
    {
        distance_from_origin = -closest_bdy_position.x();
    }

    surface_center_point = ray_direction*distance_from_origin;
    const double deltax = (max_x - min_x)/(Nxpix);//putting the points a half ... away from the edge
    const double deltay = (max_y - min_y)/(Nypix);

    ImX.resize(Nxpix*Nypix);
    ImY.resize(Nxpix*Nypix);
    I.resize(Nxpix*Nypix, geometry.parameters->nfreqs());

    threaded_for (x_idx, Nxpix,
    {
        const double xloc = min_x + deltax/2.0 + x_idx * deltax;
        for(Size y_idx = 0; y_idx<Nypix; y_idx++)
        {
              const Size totidx = y_idx + Nypix * x_idx;
              const double yloc = min_y + deltay/2.0 + y_idx * deltay;
              ImX[totidx] = xloc;
              ImY[totidx] = yloc;
        }
    });


}


//Helper function which converts surface coordinates to 3D coordinates using the saved
//Warning: use only if ImagePointPosition==ProjectionSurface
accel Vector3D Image :: surface_coords_to_3D_coordinates(const double x, const double y) const
{
    // std::cout<<"x, y: "<<x<<" "<<y<<std::endl;
    if (imagePointPosition!=ProjectionSurface)
    {
        throw std::runtime_error("Surface coordinates cannot be computed of image which does not define a projection surface.");
    }
    //TODO: maybe save these conversion things
    const double rx = ray_direction.x();
    const double ry = ray_direction.y();
    const double rz = ray_direction.z();

    const double         denominator = sqrt (rx*rx + ry*ry);
    const double inverse_denominator = 1.0 / denominator;

    const double ix =  ry * inverse_denominator;
    const double iy = -rx * inverse_denominator;

    const double jx =  rx * rz * inverse_denominator;
    const double jy =  ry * rz * inverse_denominator;
    const double jz = -denominator;

    // std::cout<<"i(x, y): "<<ix<<" "<<iy<<std::endl;
    // std::cout<<"j(x, y, z)"<<jx<<" "<<jy<<" "<<jz<<std::endl;
    // std::cout<<"denominator: "<<denominator<<std::endl;
    // std::cout<<"surface_center_point: "<<surface_center_point.x()<<" "<<surface_center_point.y()<<" "<<surface_center_point.z()<<std::endl;

    if (denominator >= 1.0e-9)
    {
        //To transform the distance in the direction back to the position, just multiply by the unit vector in that direction
        return Vector3D(x*ix + y*jx, x*iy + y*jy, y*jz) + surface_center_point;
    }
    else
    {
        return Vector3D(x, y, 0) + surface_center_point;
    }

}
