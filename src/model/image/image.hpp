#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "model/parameters/parameters.hpp"
#include "model/geometry/geometry.hpp"


enum ImageType {Intensity, OpticalDepth};


///  Image: data structure for the images
/////////////////////////////////////////
struct Image
{
    const ImageType imageType;   ///< Type of image (intensity or optical depth)
    const Size      ray_nr;      ///< number of the ray to be imaged
    const Vector3D  raydir;      ///< ray direction, if not using the

    Double1 ImX;                 ///< x coordinate of point in image
    Double1 ImY;                 ///< y coordinate of point in image
    Vector<Vector3D> ray_origin  ///< Coordinates of the start of the ray (when not using the points)

    Matrix<Real> I;              ///< intensity out along ray (index(p,f))

    Image (const Geometry& geometry, const ImageType, const Size ray_nr);
    Image (const Image& image);

    // void write (const Io &io) const;

    void set_coordinates (const Geometry& geometry);
};
