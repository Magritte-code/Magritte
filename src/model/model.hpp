#pragma once

#include "chemistry/chemistry.hpp"
#include "geometry/geometry.hpp"
#include "image/image.hpp"
#include "io/io.hpp"
#include "io/python/io_python.hpp"
#include "lines/lines.hpp"
#include "parameters/parameters.hpp"
#include "radiation/radiation.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "tools/timer.hpp"
#include "tools/types.hpp"

#include <tuple>

enum NgAccelerationType { Default, Adaptive };

struct Model {
    std::shared_ptr<Parameters> parameters;

    Geometry geometry;
    Chemistry chemistry;
    Thermodynamics thermodynamics;
    Lines lines;
    Radiation radiation;
    vector<Image> images;

    enum SpectralDiscretisation { SD_None, SD_Lines, SD_Image } spectralDiscretisation = SD_None;

    Model() :
        parameters(new Parameters()), geometry(parameters), chemistry(parameters),
        thermodynamics(parameters), lines(parameters), radiation(parameters){};

    Model(const string name) :
        parameters(new Parameters()), geometry(parameters), chemistry(parameters),
        thermodynamics(parameters), lines(parameters), radiation(parameters) {
        parameters->set_model_name(name);
        read();
    }

    void read(const Io& io);
    void write(const Io& io) const;

    void read() { read(IoPython("hdf5", parameters->model_name())); };
    void write() const { write(IoPython("hdf5", parameters->model_name())); };

    void read(const string model_name) { read(IoPython("hdf5", model_name)); };
    void write(const string model_name) const { write(IoPython("hdf5", model_name)); };

    int compute_inverse_line_widths();
    int compute_spectral_discretisation();
    int compute_spectral_discretisation(const Real width);
    int compute_spectral_discretisation(const Real nu_min, const Real nu_max);
    int compute_spectral_discretisation(
        const Real nu_min, const Real nu_max, const Size n_image_freqs);
    int compute_LTE_level_populations();
    int compute_radiation_field();
    int compute_radiation_field_feautrier_order_2();
    int compute_radiation_field_shortchar_order_0();
    int compute_Jeff();
    int compute_Jeff_sparse();
    int compute_level_populations_from_stateq();
    int compute_level_populations(const bool use_Ng_acceleration, const long max_niterations);
    int compute_level_populations_sparse(
        const bool use_Ng_acceleration, const long max_niterations);
    int compute_level_populations_shortchar(
        const bool use_Ng_acceleration, const long max_niterations);
    template <NgAccelerationType type>
    std::tuple<bool, Size> ng_acceleration_criterion(
        bool use_Ng_acceleration, Size prior_normal_iterations);
    int compute_image(const Size ray_nr);
    int compute_image_optical_depth(const Size ray_nr);
    int compute_image_new(const Vector3D raydir, const Size Nxpix,
        const Size Nypix); // actual function for the new imager
    // convenient wrappers for the new imager
    int compute_image_new(
        const double rx, const double ry, const double rz, const Size Nxpix, const Size Nypix);
    int compute_image_new(const Size ray_nr, const Size Nxpix, const Size Nypix);
    int compute_image_new(const Size ray_nr); // most similar function formulation
                                              // to old imager
    int compute_image_optical_depth_new(const Vector3D raydir, const Size Nxpix,
        const Size Nypix); // actual function for the new imager
    // convenient wrappers for the new imager
    int compute_image_optical_depth_new(
        const double rx, const double ry, const double rz, const Size Nxpix, const Size Nypix);
    int compute_image_optical_depth_new(const Size ray_nr, const Size Nxpix, const Size Nypix);
    int compute_image_optical_depth_new(const Size ray_nr); // most similar function formulation
                                                            // to old imager

    Double1 error_max;
    Double1 error_mean;

    Matrix<Real> eta;
    Matrix<Real> chi;

    Matrix<Real> S_ray;
    Matrix<Real> dtau_ray;
    Matrix<Real> u_ray;

    Matrix<Real> boundary_condition;

    Vector<Real> dshift_max;

    int set_dshift_max();

    int compute_image_for_point(const Size ray_nr, const Size p);

    int compute_radiation_field_feautrier_order_2_uv();
    int compute_radiation_field_feautrier_order_2_anis();
    int compute_radiation_field_feautrier_order_2_sparse();

    int set_eta_and_chi(const Size rr);
    int set_boundary_condition();

    Vector<Real> density;
    Matrix<Real> column;

    int set_column();
};
