#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "model/parameters/parameters.hpp"
#include "model/geometry/geometry.hpp"


enum ImageType {Intensity, OpticalDepth, PolarizedIntensity};


///  Image: data structure for the images
/////////////////////////////////////////
struct Image
{
    const ImageType imageType;   ///< Type of image (intensity or optical depth)
    const Size      ray_nr;      ///< number of the ray to be imaged

    Double1 ImX;                 ///< x coordinate of point in image
    Double1 ImY;                 ///< y coordinate of point in image

    Matrix<Real> I;              ///< intensity out along ray (index(p,f))

    Matrix<Real> I_p;            ///< intensity parallel (index(p,f))
    Matrix<Real> I_o;            ///< intensity orthogonal (index(p,f))
    Matrix<Real> U;              ///< intensity orthogonal (index(p,f))

    Image (const Geometry& geometry, const ImageType, const Size ray_nr);
    Image (const Image& image);

    // void write (const Io &io) const;

    void set_coordinates (const Geometry& geometry);
};
