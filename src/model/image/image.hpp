#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "model/parameters/parameters.hpp"
#include "model/geometry/geometry.hpp"


enum ImageType {Intensity, OpticalDepth};
enum ImagePointPosition {AllModelPoints, ProjectionSurface};
// enum ImageRayDirection

///  Image: data structure for the images
/////////////////////////////////////////
struct Image
{
    const ImageType imageType;   ///< Type of image (intensity or optical depth)
    const ImagePointPosition imagePointPosition; ///< From which points the ray starts (model points or surface outside model)
    const Size      ray_nr;      ///< number of the ray to be imaged
    const Vector3D  ray_direction;///< ray direction, if not imaging all model points
    Size closest_bdy_point;///< position of closest bdy point to start tracing ray from, if imaging a projection surface
    Vector3D surface_center_point;///< position of 0 point of projection surfa, if imaging a projection surface
    // const Vector3D  closest_bdy_point_position;///<closest bdy point to start tracing ray from, if imaging a projection surface
    // const double distance_from_origin_in_raydir;///<distance from closest boundary point to origin. Can be used to determine the 3D coordinates of the projection surface

    const Vector<Real> frequencies_to_image; ///< frequencies which need to be imaged; TODO: add this functionality

    // const bool image_all_points;///

    Double1 ImX;                 ///< x coordinate of point in image
    Double1 ImY;                 ///< y coordinate of point in image
    // Size Npixels;
    // Vector<Vector3D> ray_origin; ///< Coordinates of the start of the ray (when using the ProjectionSurface)

    Matrix<Real> I;              ///< intensity out along ray (index(p,f))


    Image (const Geometry& geometry, const ImageType it, const Size ray_nr);
    Image (const Geometry& geometry, const ImageType it, const Size ray_nr, const Size Nxpix, const Size Nypix);
    // Image (const Geometry& geometry, const ImageType it, const ImagePointPosition ipp, const Size ray_nr, const Size Nxpix, const Size Nypix);
    // template <const ImagePointPosition ipp, const Size Nxpix, const Size Nypix>
    Image (const Geometry& geometry, const ImageType it, const Vector3D ray_dir);
    Image (const Geometry& geometry, const ImageType it, const Vector3D ray_dir, const Size Nxpix, const Size Nypix);
    // Image (const Geometry& geometry, const ImageType, const Vector3D ray_dir, const bool image_all_points, const Size Nxpix, const Size Nypix);
    Image (const Image& image);

    // void write (const Io &io) const;
    // template <ImagePointPosition ipp>
    inline void set_coordinates_all_model_points (const Geometry& geometry);
    inline void set_coordinates_projection_surface (const Geometry& geometry, const Size Nxpix, const Size Nypix);

    accel Vector3D surface_coords_to_3D_coordinates(const double x, const double y) const;
    // void set_coordinates_projection_surface (const Geometry& geometry, const Size Nxpix, const Size Nypix);
};
