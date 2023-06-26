#pragma once

#include "io/io.hpp"
#include "model/geometry/geometry.hpp"
#include "model/parameters/parameters.hpp"
#include "model/radiation/frequencies/frequencies.hpp"
#include "tools/types.hpp"

enum ImageType {
    Intensity,
    OpticalDepth
}; ///< What the image represents (either intensity or optical depth)
enum ImagePointPosition {
    AllModelPoints,
    ProjectionSurface
}; ///< How to trace the imaging rays (either throughout all point or starting
   ///< from a surface outside the model)

///  Image: data structure for the images
/////////////////////////////////////////
struct Image {
    const ImageType imageType;                   ///< Type of image (intensity or optical depth)
    const ImagePointPosition imagePointPosition; ///< From which points the ray starts (model points
                                                 ///< or surface outside model)
    const Size ray_nr;                           ///< number of the ray to be imaged
    const Vector3D ray_direction;                ///< ray direction, if not imaging all model points

    // Stuff which is practically const, but the compiler complains about because
    // they are defined in a subfunction (not the constructor itself).
    Size nfreqs;                   ///< number of frequencies imaged
    Vector<Real> freqs;            ///< frequencies imaged
    Size closest_bdy_point;        ///< position of closest bdy point to start tracing
                                   ///< ray from, if imaging a projection surface
    Vector3D surface_center_point; ///< position of 0 point of projection surface,
                                   ///< if imaging a projection surface

    // The images themselves define some coordinate system, this is copied into
    // these Vector3D's
    Vector3D image_direction_x; ///< coordinate of image x direction in a 3D
                                ///< general geometry
    Vector3D image_direction_y; ///< coordinate of image y direction in a 3D
                                ///< general geometry
    Vector3D image_direction_z; ///< direction in which the image is taken

    Double1 ImX; ///< x coordinate of point in image
    Double1 ImY; ///< y coordinate of point in image

    // The actual image contents
    Matrix<Real> I; ///< value of the pixel at the corresponding coordinate

    Image(const Geometry& geometry, const Frequencies& frequencies, const ImageType it,
        const Size ray_nr);
    Image(const Geometry& geometry, const Frequencies& frequencies, const ImageType it,
        const Size ray_nr, const Size Nxpix, const Size Nypix);
    Image(const Geometry& geometry, const Frequencies& frequencies, const ImageType it,
        const Vector3D ray_dir);
    Image(const Geometry& geometry, const Frequencies& frequencies, const ImageType it,
        const Vector3D ray_dir, const Size Nxpix, const Size Nypix);
    Image(const Image& image);

    inline void set_freqs(const Frequencies& frequencies);
    void set_freqs(const Vector<Real>& image_freqs);
    inline void set_coordinates_all_model_points(const Geometry& geometry);
    inline void set_coordinates_projection_surface(
        const Geometry& geometry, const Size Nxpix, const Size Nypix);

    accel Vector3D surface_coords_to_3D_coordinates(const double x, const double y) const;
};
