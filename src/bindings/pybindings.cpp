#include "../configure.hpp" // ../ is required!
#include "io/cpp/io_cpp_text.hpp"
#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "model/parameters/parameters.hpp"
#include "pybind11/eigen.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"
#include "solver/solver.hpp"
#include "tools/types.hpp"

#include <iterator>
namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(vector<LineProducingSpecies>);
PYBIND11_MAKE_OPAQUE(vector<CollisionPartner>);

PYBIND11_MODULE(core, module) {
    // Module docstring
    module.doc() = "Core module of Magritte: a modern software library for 3D "
                   "radiative transfer.";

    // Paracabs
    module.def("pcmt_n_threads_avail", &paracabs::multi_threading::n_threads_avail,
        "Get the number of available threads (using OpenMP).");
    module.def("pcmt_set_n_threads_avail", &paracabs::multi_threading::set_n_threads_avail,
        "Set the number of available threads (using OpenMP).");
    module.def("pcmp_comm_rank", &paracabs::message_passing::comm_rank,
        "Get the rank of the current process (using MPI).");
    module.def("pcmp_comm_size", &paracabs::message_passing::comm_size,
        "Get the size of the current communicator (using MPI).");
    module.def("pcmp_length", &paracabs::message_passing::length,
        "Get the size of an array with given length that is given to the "
        "current process (using MPI).");
    module.def("pcmp_start", &paracabs::message_passing::start,
        "Get the first index of an array with given length that will be "
        "given to this process (using MPI).");
    module.def("pcmp_stop", &paracabs::message_passing::stop,
        "Get the last index of an array with given length that will be "
        "given to this process (using MPI).");

    // Define vector types
    py::bind_vector<vector<LineProducingSpecies>>(module, "vLineProducingSpecies");
    py::bind_vector<vector<CollisionPartner>>(module, "vCollisionPartner");

    // Constants
    module.attr("CC")    = CC;
    module.attr("HH")    = HH;
    module.attr("KB")    = KB;
    module.attr("AMU")   = AMU;
    module.attr("T_CMB") = T_CMB;
    // Version
    module.attr("__version__") = VERSION_NUMBER;

    // Io, base class
    py::class_<Io>(module, "Io", "Abstract input/output base class.");

    // IoText
    py::class_<IoText, Io>(
        module, "IoText", "Intput/output class for text-based io, implemented in C++.")
        // attributes
        .def_readonly(
            "io_file", &IoText::io_file, "Name of the folder or file to read from or write to.")
        // constructor
        .def(py::init<const string&>());

#if (PYTHON_IO)
    // IoPython
    py::class_<IoPython, Io>(
        module, "IoPython", "Input/output class for io, implemented in Python.")
        // attributes
        .def_readonly("implementation", &IoPython::implementation,
            "Type of input/output. Either \"hdf5\" or \"text\".")
        .def_readonly(
            "io_file", &IoPython::io_file, "Name of the folder or file to read from or write to.")
        // constructor
        .def(py::init<const string&, const string&>());
#endif

    // Solver
    py::class_<Solver>(module, "Solver")
        // constructor
        .def(py::init<>());

    // ImageType
    py::enum_<ImageType>(module, "ImageType")
        .value("Intensity", Intensity)
        .value("OpticalDepth", OpticalDepth)
        .export_values();

    // ImagePointPosition
    py::enum_<ImagePointPosition>(module, "ImagePointPosition")
        .value("AllModelPoints", AllModelPoints)
        .value("ProjectionSurface", ProjectionSurface)
        .export_values();

    // Image
    py::class_<Image>(
        module, "Image", "Image class, 2D point cloud of intensities for each frequency bin.")
        // attributes
        .def_readonly("imageType", &Image::imageType, "Type of image (intensity or optical depth).")
        .def_readonly("imagePointPosition", &Image::imagePointPosition,
            "Position of image points (model points or projection surface).")
        .def_readonly("nfreqs", &Image::nfreqs, "Number of frequency bins in the image.")
        .def_readonly("freqs", &Image::freqs, "Frequency bins of the image.")
        .def_readonly("ray_nr", &Image::ray_nr, "Number of the ray along which the image is taken.")
        .def_readonly("ImX", &Image::ImX, "X-coordinates of the points in the image plane.")
        .def_readonly("ImY", &Image::ImY, "Y-coordinates of the points in the image plane.")
        .def_readonly("image_direction_x", &Image::image_direction_x,
            "coordinate of image x direction in a 3D general geometry.")
        .def_readonly("image_direction_y", &Image::image_direction_y,
            "coordinate of image y direction in a 3D general geometry.")
        .def_readonly("image_direction_z", &Image::image_direction_z,
            "direction in which the image is taken.")
        .def_readonly(
            "I", &Image::I, "Intensity of the points in the image (for each frequency bin).")
        // constructor
        .def(py::init<const Geometry&, const Frequencies&, const ImageType&, const Size&>());

    // Model
    py::class_<Model>(module, "Model", "Class containing the Magritte model.")
        // attributes
        .def_readonly("parameters", &Model::parameters)
        .def_readonly("geometry", &Model::geometry)
        .def_readonly("chemistry", &Model::chemistry)
        .def_readonly("lines", &Model::lines)
        .def_readonly("thermodynamics", &Model::thermodynamics)
        .def_readonly("radiation", &Model::radiation)
        .def_readonly("images", &Model::images)
        .def_readwrite("eta", &Model::eta)
        .def_readwrite("chi", &Model::chi)
        .def_readonly("S_ray", &Model::S_ray)
        .def_readonly("dtau_ray", &Model::dtau_ray)
        .def_readonly("u_ray", &Model::u_ray)
        .def_readwrite("boundary_condition", &Model::boundary_condition)
        .def_readwrite("column", &Model::column)
        .def_readwrite("density", &Model::density)
        .def_readonly(
            "error_mean", &Model::error_mean, "Mean relative error in the level populations.")
        .def_readonly(
            "error_max", &Model::error_max, "Max relative error in the level populaitons.")
        // functions
        .def("read", (void(Model::*)(void)) & Model::read,
            "Read model file (using the hdf5 implementation of IoPython and "
            "using the model name specified in parameters.)")
        .def("write", (void(Model::*)(void) const) & Model::write,
            "Write model file (using the hdf5 implementation of IoPython and "
            "using the model name specified in parameters.)")
        .def("read", (void(Model::*)(const Io&)) & Model::read,
            "Read model file using the given Io object.")
        .def("write", (void(Model::*)(const Io&) const) & Model::write,
            "Write model file using the given Io object.")
        .def("read", (void(Model::*)(const string)) & Model::read,
            "Read model file (assuming HDF5 file format).")
        .def("write", (void(Model::*)(const string) const) & Model::write,
            "Write model file (assuming HDF5 file format).")
        .def("compute_inverse_line_widths", &Model::compute_inverse_line_widths,
            "Compute the inverse line widths for the model. (Needs to be "
            "recomputed whenever the temperature of turbulence changes.)")
        .def("compute_spectral_discretisation",
            (int(Model::*)(void)) & Model::compute_spectral_discretisation,
            "Compute the spectral discretisation for the model tailored for "
            "line (Gauss-Hermite) quadrature.")
        .def("compute_spectral_discretisation",
            (int(Model::*)(const Real width)) & Model::compute_spectral_discretisation,
            "Compute the spectral discretisation for the model tailored for "
            "images with the given spectral width.")
        .def("compute_spectral_discretisation",
            (int(Model::*)(const Real nu_min, const Real nu_max))
                & Model::compute_spectral_discretisation,
            "Compute the spectral discretisation for the model tailored for "
            "images with the given min and max frequency.")
        .def("compute_spectral_discretisation",
            (int(Model::*)(const Real nu_min, const Real nu_max, const Size n_image_freqs))
                & Model::compute_spectral_discretisation,
            "Compute the spectral discretisation for the model tailored for "
            "images with the given min and max frequency. Can also specify the "
            "amount of frequency bins to use (instead of defaulting to "
            "parameters.nfreqs).")
        .def("compute_LTE_level_populations", &Model::compute_LTE_level_populations,
            "Compute the level populations for the model assuming local "
            "thermodynamic equilibrium (LTE).")
        .def("compute_radiation_field_feautrier_order_2",
            &Model::compute_radiation_field_feautrier_order_2,
            "Compute the radiation field for the modle using the 2nd-order "
            "Feautrier solver.")
        .def("compute_radiation_field_feautrier_order_2_sparse",
            &Model::compute_radiation_field_feautrier_order_2_sparse,
            "Compute the radiation field for the modle using the 2nd-order "
            "Feautrier solver.")
        /// Solver is bugged, so removed from the api, as the shortchar solver can
        /// replace it
        .def("compute_radiation_field_feautrier_order_2_uv",
            &Model::compute_radiation_field_feautrier_order_2_uv,
            "Compute the radiation field for the modle using the 2nd-order Feautrier solver.")
        .def("compute_radiation_field_feautrier_order_2_anis",
            &Model::compute_radiation_field_feautrier_order_2_anis,
            "Compute the radiation field for the modle using the 2nd-order "
            "Feautrier solver, anisotropic case.")
        .def("compute_radiation_field_shortchar_order_0",
            &Model::compute_radiation_field_shortchar_order_0,
            "Compute the radiation field for the modle using the 0th-order "
            "short-characteristics methods.")
        .def("compute_Jeff", &Model::compute_Jeff,
            "Compute the effective mean intensity in the line.")
        .def("compute_level_populations_from_stateq", &Model::compute_level_populations_from_stateq,
            "Compute the level populations for the model assuming statistical "
            "equilibrium.")
        .def("compute_level_populations", &Model::compute_level_populations,
            "Compute the level populations for the model assuming statistical "
            "equilibrium until convergence, optionally using Ng-acceleration, "
            "and for the given maximum number of iterations.")
        .def("compute_level_populations_sparse", &Model::compute_level_populations_sparse,
            "Compute the level populations for the model assuming statistical "
            "equilibrium until convergence, optionally using Ng-acceleration, "
            "and for the given maximum number of iterations. (Memory sparse "
            "option.)")
        .def("compute_level_populations_shortchar", &Model::compute_level_populations_shortchar,
            "Compute the level populations using the short-characteristics "
            "solver for the model assuming statistical equilibrium until "
            "convergence, optionally using Ng-acceleration, and for the given "
            "maximum number of iterations.")
        .def("compute_image", &Model::compute_image,
            "Compute an image for the model along the given ray.")
        .def("compute_image_new", (int(Model::*)(const Size ray_nr)) & Model::compute_image_new,
            "Compute an image of the model along the given ray direction, using "
            "the new imager.")
        .def("compute_image_new",
            (int(Model::*)(const Size ray_nr, const Size Nxpix, const Size Nypix))
                & Model::compute_image_new,
            "Compute an image of the model along the given ray direction, using "
            "the new imager, specifying the image resolution.")
        .def("compute_image_new",
            (int(Model::*)(const double rx, const double ry, const double rz, const Size Nxpix,
                const Size Nypix))
                & Model::compute_image_new,
            "Compute an image of the model along the given ray direction, using "
            "the new imager, specifying the ray direction and image resolution.")
        .def("compute_image_optical_depth", &Model::compute_image_optical_depth,
            "Compute an image of the optical depth for the model along the "
            "given ray.")
        .def("compute_image_optical_depth_new",
            (int(Model::*)(const Size ray_nr)) & Model::compute_image_optical_depth_new,
            "Compute an image of the optical depth for the model along the "
            "given ray direction, using the new imager.")
        .def("compute_image_optical_depth_new",
            (int(Model::*)(const Size ray_nr, const Size Nxpix, const Size Nypix))
                & Model::compute_image_optical_depth_new,
            "Compute an image of the optical depth for the model along the "
            "given ray direction, using the new imager, specifying the image "
            "resolution.")
        .def("compute_image_optical_depth_new",
            (int(Model::*)(const double rx, const double ry, const double rz, const Size Nxpix,
                const Size Nypix))
                & Model::compute_image_optical_depth_new,
            "Compute an image of the optical depth for the model along the "
            "given ray direction, using the new imager, specifying the ray "
            "direction and image resolution.")
        .def("set_eta_and_chi", &Model::set_eta_and_chi,
            "Set latest emissivity and opacity for the model in the eta and chi "
            "variables respectively.")
        .def("set_boundary_condition", &Model::set_boundary_condition,
            "Set boundary condition (internally).")
        .def("compute_image_for_point", &Model::compute_image_for_point,
            "Compute image (single pixel) for a single point.")
        .def("set_column", &Model::set_column, "Set column (internally).")
        // constructor
        .def(py::init<const string>())
        .def(py::init<>());

    // Parameters
    py::class_<Parameters, std::shared_ptr<Parameters>>(
        module, "Parameters", "Class containing the model parameters.")
        // io
        .def_readwrite("n_off_diag", &Parameters::n_off_diag,
            "Bandwidth of the ALO (0=diagonal, 1=tri-diaginal, "
            "2=penta-diagonal, ...)")
        .def_readwrite("max_width_fraction", &Parameters::max_width_fraction,
            "Max tolerated Doppler shift as fraction of line width "
            "(default=0.5).")
        .def_readwrite("convergence_fraction", &Parameters::convergence_fraction,
            "Fraction of levels that should obey the convergence criterion.")
        .def_readwrite("min_rel_pop_for_convergence", &Parameters::min_rel_pop_for_convergence,
            "Minimum relative level population to be considered in "
            "the convergence criterion.")
        .def_readwrite("pop_prec", &Parameters::pop_prec, "Required precision for ALI.")
        .def_readwrite("min_opacity", &Parameters::min_opacity,
            "Minimum opacity that will be assumed in the solver.")
        .def_readwrite("min_line_opacity", &Parameters::min_line_opacity,
            "Minimum line opacity that will be assumed in the solver.")
        .def_readwrite("min_dtau", &Parameters::min_dtau,
            "Minimum optical depth increment that will be assumed in the solver.")
        .def_readwrite("population_inversion_fraction", &Parameters::population_inversion_fraction,
            "Threshold factor for population inversion required for LTE to be used; set this "
            "higher than 1")
        .def_readwrite("store_intensities", &Parameters::store_intensities,
            "Whether or not to store intensities.")
        .def_readwrite("one_line_approximation", &Parameters::one_line_approximation,
            "Whether or not to use one line approximation.")
        .def_readwrite("sum_opacity_emissivity_over_all_lines",
            &Parameters::sum_opacity_emissivity_over_all_lines,
            "Whether or not to sum the opacity, emissivity over all "
            "lines (including the zero contributions).")
        .def_readwrite("max_distance_opacity_contribution",
            &Parameters::max_distance_opacity_contribution,
            "The distance (scaled to line width) at which we ignore "
            "the opacity/emissivity contribution.")
        .def_readwrite("use_adaptive_Ng_acceleration", &Parameters::use_adaptive_Ng_acceleration,
            "Whether to use adaptive Ng acceleration.")
        .def_readwrite("Ng_acceleration_mem_limit", &Parameters::Ng_acceleration_mem_limit,
            "Memory limit for ng acceleration.")
        .def_readwrite("adaptive_Ng_acceleration_use_max_criterion",
            &Parameters::adaptive_Ng_acceleration_use_max_criterion,
            "Whether of not to use the maximum relative change as "
            "criterion for the adaptive ng-acceleration")
        .def_readwrite("Ng_acceleration_remove_N_its", &Parameters::Ng_acceleration_remove_N_its,
            "Number of iterations to ignore when using ng-acceleration")
        .def_readwrite("adaptive_Ng_acceleration_min_order",
            &Parameters::adaptive_Ng_acceleration_min_order,
            "Minimal order of ng acceleration when using adaptive Ng "
            "acceleration. Has to be larger than 1.")
        // setters
        .def("set_model_name", &Parameters::set_model_name, "Set model name.")
        .def("set_dimension", &Parameters::set_dimension, "Set spatial dimension of the model.")
        .def("set_npoints", &Parameters::set_npoints, "Set number of points.")
        .def("set_nrays", &Parameters::set_nrays, "Set number of rays.")
        .def("set_hnrays", &Parameters::set_hnrays, "Set the half of the number of rays.")
        .def("set_nboundary", &Parameters::set_nboundary, "Set number of boundary points.")
        .def("set_nfreqs", &Parameters::set_nfreqs,
            "Set number of frequency bins.") // TODO: phase this out, as this is
                                             // a) not essential to model setup
                                             // and b) dangerous to play around
                                             // with (segfaults if wrong value is
                                             // set)
        .def("set_nspecs", &Parameters::set_nspecs, "Set number of species.")
        .def("set_nlspecs", &Parameters::set_nlspecs, "Set number of line producing species.")
        .def("set_nlines", &Parameters::set_nlines, "Set number of lines.")
        .def("set_nquads", &Parameters::set_nquads, "Set number of quadrature points.")
        .def("set_use_scattering", &Parameters::set_use_scattering,
            "Set whether or not to use scattering.")
        .def("set_spherical_symmetry", &Parameters::set_spherical_symmetry,
            "Set whether or not to use spherical symmetry")
        .def("set_adaptive_ray_tracing", &Parameters::set_adaptive_ray_tracing,
            "Set whether or not to use adaptive ray tracing.")
        // getters
        .def("version", &Parameters::version, "Magritte version number.")
        .def("model_name", &Parameters::model_name, "Model name.")
        .def("dimension", &Parameters::dimension, "Spatial dimension of the model.")
        .def("npoints", &Parameters::npoints, "Number of points.")
        .def("nrays", &Parameters::nrays, "Number of rays.")
        .def("hnrays", &Parameters::hnrays, "Half the number of rays.")
        .def("nboundary", &Parameters::nboundary, "Number of boundary points.")
        .def("nfreqs", &Parameters::nfreqs, "Number of frequency bins.")
        .def("nspecs", &Parameters::nspecs, "Number of species.")
        .def("nlspecs", &Parameters::nlspecs, "Number of line producing species.")
        .def("nlines", &Parameters::nlines, "Number of lines.")
        .def("nquads", &Parameters::nquads, "Number of quadrature points.")
        .def("use_scattering", &Parameters::use_scattering, "Whether or not to use scattering.")
        .def("spherical_symmetry", &Parameters::spherical_symmetry,
            "Whether or not to use spherical symmetry.")
        .def("adaptive_ray_tracing", &Parameters::adaptive_ray_tracing,
            "Whether or not to use adaptive ray tracing.")
        // functions
        .def("read", &Parameters::read, "Rread object from file.")
        .def("write", &Parameters::write, "Write object to file.");

    // Geometry
    py::class_<Geometry>(module, "Geometry", "Class containing the model geometry.")
        // attributes
        .def_readonly("points", &Geometry::points, "Points object.")
        .def_readonly("rays", &Geometry::rays, "Rays object.")
        .def_readonly("boundary", &Geometry::boundary, "Boundary object")
        .def_readonly("lengths", &Geometry::lengths,
            "Array containing the lengths of the rays for each "
            "direction and point.")
        // io
        .def("read", &Geometry::read, "Read object from file.")
        .def("write", &Geometry::write, "Write object to file.");

    // Points
    py::class_<Points>(module, "Points", "Class containing the spatial points.")
        // attributes
        .def_readwrite("position", &Points::position, "Array with position vectors of the points.")
        .def_readwrite("velocity", &Points::velocity,
            "Array with velocity vectors of the points (as a fraction "
            "of the speed of light).")
        .def_readwrite(
            "neighbors", &Points::neighbors, "Linearised array of neighbours of each point.")
        .def_readwrite("n_neighbors", &Points::n_neighbors, "Number of neighbours of each point.")
        .def_readonly("cum_n_neighbors", &Points::cum_n_neighbors,
            "Cumulative number of neighbours of each point.")
        // io
        .def("read", &Points::read, "Read object from file.")
        .def("write", &Points::write, "Write object to file.");

    // Rays
    py::class_<Rays>(
        module, "Rays", "Class containing the (light) rays (directional discretisation).")
        // attributes
        .def_readwrite("direction", &Rays::direction, "Array with direction vector of each ray.")
        .def_readonly(
            "antipod", &Rays::antipod, "Array with the number of the antipodal ray for each ray.")
        .def_readwrite("weight", &Rays::weight,
            "Array with the weights that each ray contributes in "
            "integrals over directions.")
        // io
        .def("read", &Rays::read, "Read object from file.")
        .def("write", &Rays::write, "Write object to file.");

    // Boundary Condition
    py::enum_<BoundaryCondition>(module, "BoundaryCondition")
        .value("Zero", Zero)
        .value("Thermal", Thermal)
        .value("CMB", CMB)
        .export_values();

    // Boundary
    py::class_<Boundary>(module, "Boundary", "Class containing the model boundary.")
        // attributes
        .def_readwrite("boundary2point", &Boundary::boundary2point,
            "Array with point index for each boundary point.")
        .def_readonly("point2boundary", &Boundary::point2boundary,
            "Array with boundary index for each point.")
        .def_readwrite("boundary_temperature", &Boundary::boundary_temperature,
            "Array with radiative temperature for each boundary point "
            "(only relevant for thermal boundary conditions).")
        // functions
        .def("set_boundary_condition", &Boundary::set_boundary_condition,
            "Setter for the boundary condition.")
        .def("get_boundary_condition", &Boundary::get_boundary_condition,
            "Getter for the boundary condition.")
        // io
        .def("read", &Boundary::read, "Read object from file.")
        .def("write", &Boundary::write, "Write object to file.");

    // Thermodynamics
    py::class_<Thermodynamics>(module, "Thermodynamics", "Class containing the thermodynamics.")
        // attributes
        .def_readonly("temperature", &Thermodynamics::temperature, "Temperature object.")
        .def_readonly("turbulence", &Thermodynamics::turbulence, "Turbulence object.")
        // io
        .def("read", &Thermodynamics::read, "Read object from file.")
        .def("write", &Thermodynamics::write, "Write object to file.");

    // Temperature
    py::class_<Temperature>(module, "Temperature", "Class containing the temperature.")
        // attributes
        .def_readwrite("gas", &Temperature::gas, "Kinetic temperature of the gas.")
        // functions
        .def("read", &Temperature::read, "Read object from file.")
        .def("write", &Temperature::write, "Write object to file.");

    // Turbulence
    py::class_<Turbulence>(module, "Turbulence", "Class containing the (micro) turbulence.")
        // attributes
        .def_readwrite("vturb2", &Turbulence::vturb2,
            "Square of the micro turbulence as a fraction of the speed of light.")
        // functions
        .def("read", &Turbulence::read, "Read object from file.")
        .def("write", &Turbulence::write, "Write object to file.");

    // Chemistry
    py::class_<Chemistry>(module, "Chemistry", "Class containing the chemistry.")
        // attributes
        .def_readonly("species", &Chemistry::species, "Species object.")
        // functions
        .def("read", &Chemistry::read, "Read object from file.")
        .def("write", &Chemistry::write, "Write object to file.");

    // Species
    py::class_<Species>(module, "Species", "Class containing the chemical species.")
        // attributes
        .def_readwrite("symbol", &Species::symbol, "Symbol of the species.")
        .def_readwrite("abundance", &Species::abundance, "Array with the abundances at each point.")
        // functions
        .def("read", &Species::read, "Read object from file.")
        .def("write", &Species::write, "Write object to file.");

    // Lines
    py::class_<Lines>(module, "Lines", "Class containing the lines and their data.")
        // attributes
        .def_readonly(
            "lineProducingSpecies", &Lines::lineProducingSpecies, "Vector of LineProducingSpecies.")
        .def_readonly("emissivity", &Lines::emissivity,
            "Array with emissivities for each point and frequency bin.")
        .def_readonly(
            "opacity", &Lines::opacity, "Array with opacities for each point and frequency bin.")
        .def_readonly("inverse_width", &Lines::inverse_width,
            "Array with the inverse widths for each line at each point.")
        .def_readonly("line", &Lines::line, "Array with the line (centre) frequencies.")
        // functions
        .def("read", &Lines::read, "Read object from file.")
        .def("write", &Lines::write, "Write object to file.")
        .def("set_emissivity_and_opacity", &Lines::set_emissivity_and_opacity,
            "Set the emissivity and opacity arrays.")
        .def("resize_LineProducingSpecies", &Lines::resize_LineProducingSpecies,
            "Resize the vector LineProducingSpecies.");

    // LineProducingSpecies
    py::class_<LineProducingSpecies>(
        module, "LineProducingSpecies", "Class containing a line producing species.")
        // attributes
        .def_readonly("linedata", &LineProducingSpecies::linedata, "Linedata object.")
        .def_readonly("quadrature", &LineProducingSpecies::quadrature, "Quadrature object.")
        .def_readonly("Lambda", &LineProducingSpecies::lambda,
            "Lambda object") // "lambda" is invalid in Python, use "Lambda"
        .def_readwrite("Jeff", &LineProducingSpecies::Jeff,
            "Array with effective mean intensities in the lines.")
        .def_readwrite("Jlin", &LineProducingSpecies::Jlin)
        .def_readonly("Jdif", &LineProducingSpecies::Jdif)
        .def_readonly("fraction_not_converged", &LineProducingSpecies::fraction_not_converged)

        .def_readonly("J", &LineProducingSpecies::J, "Isotropic radiation field.")
        .def_readonly(
            "J2_0", &LineProducingSpecies::J2_0, "Anisotropic radiation field tensor element 0")
        .def_readonly("J2_1_Re", &LineProducingSpecies::J2_1_Re,
            "Anisotropic radiation field tensor element 1 (real part)")
        .def_readonly("J2_1_Im", &LineProducingSpecies::J2_1_Im,
            "Anisotropic radiation field tensor element 1 (imaginary part)")
        .def_readonly("J2_2_Re", &LineProducingSpecies::J2_2_Re,
            "Anisotropic radiation field tensor element 2 (real part)")
        .def_readonly("J2_2_Im", &LineProducingSpecies::J2_2_Im,
            "Anisotropic radiation field tensor element 2 (imaginary part)")
        .def_readwrite("nr_line", &LineProducingSpecies::nr_line)
        .def_readwrite("population", &LineProducingSpecies::population,
            "Array with level populations for each point.")
        .def_readwrite("population_tot", &LineProducingSpecies::population_tot,
            "Array with the sum of all level populations at each point. (Should "
            "be equal to the abundance of the species.)")
        .def_readonly("population_prev1", &LineProducingSpecies::population_prev1)
        .def_readonly("population_prev2", &LineProducingSpecies::population_prev2)
        .def_readonly("population_prev3", &LineProducingSpecies::population_prev3)
        .def_readonly("populations", &LineProducingSpecies::populations)
        .def_readonly("RT", &LineProducingSpecies::RT)
        .def_readonly("LambdaStar", &LineProducingSpecies::LambdaStar)
        .def_readonly("LambdaTest", &LineProducingSpecies::LambdaTest)
        // functions
        .def("read", &LineProducingSpecies::read, "Read object from file.")
        .def("write", &LineProducingSpecies::write, "Write object to file.")
        .def("index", &LineProducingSpecies::index);

    // Lambda
    py::class_<Lambda>(module, "Lambda")
        // attributes
        .def_readonly("Ls", &Lambda::Ls)
        .def_readonly("nr", &Lambda::nr)
        // .def_readwrite ("size", &Lambda::size)
        .def_readonly("Lss", &Lambda::Lss)
        .def_readonly("nrs", &Lambda::nrs)
        // functions
        .def("add_element", &Lambda::add_element)
        .def("linearize_data", &Lambda::linearize_data)
        .def("MPI_gather", &Lambda::MPI_gather);

    // Quadrature
    py::class_<Quadrature>(
        module, "Quadrature", "Class containing the data for Gauss-Hermite quadrature.")
        // attributes
        .def_readwrite("roots", &Quadrature::roots,
            "Array containing the roots for the Gauss-Hermite quadrature.")
        .def_readwrite("weights", &Quadrature::weights,
            "Array containing the weights for the Gauss-Hermite quadrature.")
        // functions
        .def("read", &Quadrature::read, "Read object from file.")
        .def("write", &Quadrature::write, "Write object to file.");

    // Linedata
    py::class_<Linedata>(module, "Linedata", "Class containing line data.")
        // attributes
        .def_readwrite("num", &Linedata::num, "Number of the species corresponding the line data.")
        .def_readwrite("sym", &Linedata::sym, "Chemical symbol of the species.")
        .def_readwrite(
            "inverse_mass", &Linedata::inverse_mass, "Inverse mass (1/mass) of the species.")
        .def_readwrite("nlev", &Linedata::nlev, "Number of energy levels.")
        .def_readwrite("nrad", &Linedata::nrad, "Number of radiative transitions.")
        .def_readwrite(
            "irad", &Linedata::irad, "Array with upper levels of the radiative transitions.")
        .def_readwrite(
            "jrad", &Linedata::jrad, "Array with lower levels of the radiative transitions.")
        .def_readwrite("energy", &Linedata::energy, "Array with the energy for each level.")
        .def_readwrite(
            "weight", &Linedata::weight, "Array with the statistical weight for each level.")
        .def_readwrite(
            "frequency", &Linedata::frequency, "Array with frequencies of the line transitions.")
        .def_readwrite("A", &Linedata::A, "Array with Einstein A coefficients.")
        .def_readwrite("Ba", &Linedata::Ba, "Array with Einstein B (absorption) coeficients.")
        .def_readwrite(
            "Bs", &Linedata::Bs, "Array with Einstsin B (stimulated emission) coefficients.")
        .def_readwrite("ncolpar", &Linedata::ncolpar, "Number of collision partners.")
        .def_readwrite("colpar", &Linedata::colpar, "Vector of collision partner objects.")
        // functions
        .def("read", &Linedata::read, "Read object from file.")
        .def("write", &Linedata::write, "Write object to file.")
        // constructor
        .def(py::init<>());

    // Colpartner
    py::class_<CollisionPartner>(
        module, "CollisionPartner", "Class containing collision partner data.")
        // attributes
        .def_readwrite("num_col_partner", &CollisionPartner::num_col_partner,
            "Number of the species corresponding to the collision partner.")
        .def_readwrite("orth_or_para_H2", &CollisionPartner::orth_or_para_H2,
            "In case the collision partner is H2, this indicates "
            "whether it is ortho (o) or para (p).")
        .def_readwrite(
            "ntmp", &CollisionPartner::ntmp, "Number of temperatures at which data is given.")
        .def_readwrite("ncol", &CollisionPartner::ncol, "Number of collisional transitions.")
        .def_readwrite("icol", &CollisionPartner::icol,
            "Array with upper levels of the collisional transitions.")
        .def_readwrite("jcol", &CollisionPartner::jcol,
            "Array with lower levels of the collisional transitions.")
        .def_readwrite("tmp", &CollisionPartner::tmp,
            "Array with temperatures corresponding to the collisional data.")
        .def_readwrite("Ce", &CollisionPartner::Ce, "Array with collisional excitation rates.")
        .def_readwrite("Cd", &CollisionPartner::Cd, "Array with collisional de-excitation rates.")
        // functions
        .def("read", &CollisionPartner::read, "Read object from file.")
        .def("write", &CollisionPartner::write, "Write object to file.")
        // constructor
        .def(py::init<>());

    // Radiation
    py::class_<Radiation>(module, "Radiation", "Class containing the radiation field.")
        // attributes
        .def_readonly("frequencies", &Radiation::frequencies, "Frequencies object.")
        .def_readonly("I", &Radiation::I,
            "Array containing the intensity for each ray, point, and "
            "frequency bin.")
        .def_readonly("u", &Radiation::u,
            "Array containing the mean intensity up and down a ray for "
            "each ray pair, point, and frequency bin.")
        .def_readonly("v", &Radiation::v,
            "Array containing the flux up and down a ray for each ray "
            "pair, point, and frequency bin.")
        .def_readonly("J", &Radiation::J,
            "Array containing the mean intensity for each point and "
            "frequency bin.")
        // functions
        .def("read", &Radiation::read, "Read object from file.")
        .def("write", &Radiation::write, "Write object to file.");

    // Frequencies
    py::class_<Frequencies>(module, "Frequencies", "Class containing the (local) frequency bins.")
        // attributes
        .def_readonly(
            "nu", &Frequencies::nu, "Array of frequency bins for each point in the model.")
        // functions
        .def("read", &Frequencies::read, "Read object from file.")
        .def("write", &Frequencies::write, "Write object to file.");

    // Vector <Size>
    py::class_<Vector<Size>>(module, "VSize", py::buffer_protocol())
        // buffer
        .def_buffer([](Vector<Size>& v) -> py::buffer_info {
            return py::buffer_info(v.vec.data(),       // Pointer to buffer
                sizeof(Size),                          // Size of one element
                py::format_descriptor<Size>::format(), // Python struct-style format
                                                       // descriptor
                1,                                     // Number of dimensions
                {v.vec.size()},                        // Buffer dimensions
                {sizeof(Size)}                         // Strides (in bytes) for each index
            );
        })
        // functions
        .def("set", &Vector<Size>::set_1D_array)
        // constructor
        .def(py::init());

    // Vector <Real>
    py::class_<Vector<Real>>(module, "VReal", py::buffer_protocol())
        // buffer
        .def_buffer([](Vector<Real>& v) -> py::buffer_info {
            return py::buffer_info(v.vec.data(),       // Pointer to buffer
                sizeof(Real),                          // Size of one element
                py::format_descriptor<Real>::format(), // Python struct-style format
                                                       // descriptor
                1,                                     // Number of dimensions
                {v.vec.size()},                        // Buffer dimensions
                {sizeof(Real)}                         // Strides (in bytes) for each index
            );
        })
        // functions
        .def("set", &Vector<Real>::set_1D_array)
        // constructor
        .def(py::init());

    // Vector3D
    py::class_<Vector3D>(module, "V3DReal", py::buffer_protocol())
        // buffer
        .def_buffer([](Vector3D& v) -> py::buffer_info {
            return py::buffer_info(&v.data[0],         // Pointer to buffer
                sizeof(Real),                          // Size of one element
                py::format_descriptor<Real>::format(), // Python struct-style format
                                                       // descriptor
                1,                                     // Number of dimensions
                {3},                                   // Buffer dimensions
                {sizeof(Real)}                         // Strides (in bytes) for each index
            );
        })
        // functions TODO IMPLEMENT when necessary (for now, no Vector3D is
        // expected* to be set by the user). *We might in the future try to ask
        // for a direction in '(numpy vector to) Vector3D' format .def ("set",
        // &Vector3D<Real>::set_1D_array) constructor
        .def(py::init());

    // Matrix <Real>
    py::class_<Matrix<Real>, Vector<Real>>(module, "MReal", py::buffer_protocol())
        // buffer
        .def_buffer([](Matrix<Real>& m) -> py::buffer_info {
            return py::buffer_info(m.vec.data(),       // Pointer to buffer
                sizeof(Real),                          // Size of one element
                py::format_descriptor<Real>::format(), // Python struct-style format
                                                       // descriptor
                2,                                     // Number of dimensions
                py::detail::any_container<ssize_t>({m.nrows, m.ncols}), // Buffer dimensions
                py::detail::any_container<ssize_t>(
                    {sizeof(Real) * m.ncols, sizeof(Real)}) // Strides (in bytes) for each index
            );
        })
        .def_readwrite("vec", &Vector<Real>::vec)
        .def_readwrite("nrows", &Matrix<Real>::nrows)
        .def_readwrite("ncols", &Matrix<Real>::ncols)
        // functions
        .def("set", &Matrix<Real>::set_2D_array)
        // constructor
        .def(py::init());

    // Matrix <Size>
    py::class_<Matrix<Size>, Vector<Size>>(module, "MSize", py::buffer_protocol())
        // buffer
        .def_buffer([](Matrix<Size>& m) -> py::buffer_info {
            return py::buffer_info(m.vec.data(),       // Pointer to buffer
                sizeof(Size),                          // Size of one element
                py::format_descriptor<Size>::format(), // Python struct-style format
                                                       // descriptor
                2,                                     // Number of dimensions
                py::detail::any_container<ssize_t>({m.nrows, m.ncols}), // Buffer dimensions
                py::detail::any_container<ssize_t>(
                    {sizeof(Size) * m.ncols, sizeof(Size)}) // Strides (in bytes) for each index
            );
        })
        .def_readwrite("vec", &Vector<Size>::vec)
        .def_readwrite("nrows", &Matrix<Size>::nrows)
        .def_readwrite("ncols", &Matrix<Size>::ncols)
        // functions
        .def("set", &Matrix<Size>::set_2D_array)
        // constructor
        .def(py::init());

    // Tensor <Real>
    py::class_<Tensor<Real>, Vector<Real>>(module, "TReal", py::buffer_protocol())
        // buffer
        .def_buffer([](Tensor<Real>& t) -> py::buffer_info {
            return py::buffer_info(t.vec.data(),       // Pointer to buffer
                sizeof(Real),                          // Size of one element
                py::format_descriptor<Real>::format(), // Python struct-style format
                                                       // descriptor
                3,                                     // Number of dimensions
                py::detail::any_container<ssize_t>(
                    {t.nrows, t.ncols, t.depth}), // Buffer dimensions
                py::detail::any_container<ssize_t>({sizeof(Real) * t.ncols * t.depth,
                    sizeof(Real) * t.depth, sizeof(Real)}) // Strides (in bytes) for each index
            );
        })
        .def_readwrite("vec", &Vector<Real>::vec)
        .def_readwrite("nrows", &Tensor<Real>::nrows)
        .def_readwrite("ncols", &Tensor<Real>::ncols)
        .def_readwrite("depth", &Tensor<Real>::depth)
        // functions
        .def("set", &Tensor<Real>::set_3D_array)
        // constructor
        .def(py::init());

    // Vector <Vector3D>
    py::class_<Vector<Vector3D>>(module, "VVector3D", py::buffer_protocol())
        // buffer
        .def_buffer([](Vector<Vector3D>& v) -> py::buffer_info {
            return py::buffer_info(v.vec.data(),         // Pointer to buffer
                sizeof(double),                          // Size of one element
                py::format_descriptor<double>::format(), // Python struct-style
                                                         // format descriptor
                2,                                       // Number of dimensions
                py::detail::any_container<ssize_t>(
                    {(ssize_t)v.vec.size(), (ssize_t)3}), // Buffer dimensions
                py::detail::any_container<ssize_t>(
                    {sizeof(double) * 3, sizeof(double)}) // Strides (in bytes) for each index
            );
        })
        // functions
        .def("set", (void(Vector<Vector3D>::*)(
                        py::array_t<double, py::array::c_style | py::array::forcecast>))
                        & Vector<Vector3D>::set_2D_array)
        // constructor
        .def(py::init());
}

// namespace pybind11
//{
//     namespace detail
//     {
//         template<> struct type_caster<Vector3D>
//         {
//             PYBIND11_TYPE_CASTER(Vector3D, _("Vector3D"));
//
//             // Conversion part 1 (Python -> C++)
//             bool load(py::handle src, bool convert)
//             {
//                 // if ( !convert and !py::array_t<T>::check_(src) )
//                     // return false;
//
//                 // auto buf = py::array_t<Real, py::array::c_style |
//                 py::array::forcecast>::ensure(src); auto buf =
//                 array::ensure(src);
//
//                 if ( !buf )
//                     return false;
//
//                 auto dims = buf.ndim();
//                 if ( dims != 1  )
//                     return false;
//
//                 // std::vector<size_t> shape(1);
//                 // shape[0] = buf.shape()[0];
//                 Real* data = (Real*) buf.data();
//
//                 value = Vector3D (data[0], data[1], data[2]);
//
//                 // auto ref =
//                 reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                 // memcpy(value.vec.data(), buf.data(),
//                 shape[0]*sizeof(Vector3D));
//
//                 // int result =
//                 detail::npy_api::get().PyArray_CopyInto_(ref.ptr(),
//                 buf.ptr());
//
//                 // if (result < 0)   // Copy failed!
//                 // {
//                    // PyErr_Clear();
//                    // return false;
//                 // }
//
//                 return true;
//             }
//
//             //Conversion part 2 (C++ -> Python)
//             static py::handle cast(const Vector3D& src,
//             py::return_value_policy policy, py::handle parent)
//             {
//                 std::vector<size_t> shape(1);
//                 shape[0] = 3;
//
//                 std::vector<size_t> strides(1);
//                 strides[0] = sizeof(Real);
//
//                 py::array arr (std::move(shape), std::move(strides),
//                 &(src.data[0]));
//
//                 return arr.release();
//             }
//         };
//
//
//         template<> struct type_caster<Vector<Vector3D>>
//         {
//             PYBIND11_TYPE_CASTER(Vector<Vector3D>, _("Vector<Vector3D>"));
//
//             // Conversion part 1 (Python -> C++)
//             bool load(py::handle src, bool convert)
//             {
//                 auto buf = array::ensure(src);
//                 if (!buf) return false;
//
//                 auto dims = buf.ndim();
//                 if (dims != 2) return false;
//
//                 std::vector<size_t> shape(2);
//                 shape[0] = buf.shape()[0];
//                 shape[1] = buf.shape()[1];
//
//                 if (shape[1] != 3) return false;
//
//                 // value = Vector<Vector3D> (shape[0]);
//                 value.resize(shape[0]);
//
//                 // Vector3D>* data = (Vector<Vector3D>*) buf.data();
//
//                 // for (size_t i = 0; i < shape[0]; i++)
//                 // {
//                     // value[i] = Vector3D()
//                 // }
//
//                 // Vector<Vector3D>* data = (Vector<Vector3D>*) buf.data();
//
//                 // auto ref =
//                 py::reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                 memcpy(&(value.vec[0].data[0]), buf.data(),
//                 shape[0]*shape[1]*sizeof(Real));
//
//                 // int result =
//                 detail::npy_api::get().PyArray_CopyInto_(ref.ptr(),
//                 buf.ptr());
//
//                 // if (result < 0)   // Copy failed!
//                 // {
//                    // PyErr_Clear();
//                    // return false;
//                 // }
//
//                 return true;
//             }
//
//             //Conversion part 2 (C++ -> Python)
//             static py::handle cast(const Vector<Vector3D>& src,
//             py::return_value_policy policy, py::handle parent)
//             {
//                 std::vector<size_t> shape(2);
//                 shape[0] = src.vec.size();
//                 shape[1] = 3;
//
//                 std::vector<size_t> strides(2);
//                 strides[0] = sizeof(Real)*shape[1];
//                 strides[1] = sizeof(Real);
//
//                 py::array arr (std::move(shape), std::move(strides),
//                 &(src.vec[0].data[0]));
//
//                 return arr.release();
//             }
//         };
//
//
//         template<typename type>
//         struct type_caster<Vector<type>>
//         {
//             PYBIND11_TYPE_CASTER(Vector<type>, _("Vector<type>"));
//
//             // Conversion part 1 (Python -> C++)
//             bool load(py::handle src, bool convert)
//             {
//                 // if ( !convert and !py::array_t<T>::check_(src) )
//                     // return false;
//
//                 // auto buf = py::array_t<Real, py::array::c_style |
//                 py::array::forcecast>::ensure(src); auto buf =
//                 py::array_t<type, py::array::c_style |
//                 py::array::forcecast>::ensure(src);
//
//                 if ( !buf )
//                     return false;
//
//                 auto dims = buf.ndim();
//                 if ( dims != 1  )
//                     return false;
//
//                 std::vector<size_t> shape(1);
//                 shape[0] = buf.shape()[0];
//
//                 cout << "----------- value.data() IN  " << value.vec.data()
//                 << endl; cout << "----------- value.dat    IN  " << value.dat
//                 << endl;
//
//                 value = Vector<type> (shape[0]);
//
//                 // std::copy (buf.data(), buf.data()+buf.size(),
//                 value.vec.data()); cout << "buf.size() = " << buf.size() <<
//                 endl; cout << "shape[0]   = " << shape[0]   << endl; cout <<
//                 "----------- value.data() IN  " << value.vec.data() << endl;
//                 cout << "----------- value.dat    IN  " << value.dat << endl;
//                 value.set_dat();
//                 cout << "----------- value.data() IN  " << value.vec.data()
//                 << endl; cout << "----------- value.dat    IN  " << value.dat
//                 << endl; std::copy (buf.data(), buf.data()+buf.size(),
//                 value.dat); cout << "----------- value.data() OUT " <<
//                 value.vec.data() << endl; cout << "----------- value.dat OUT
//                 " << value.dat        << endl;
//
//                 // auto ref =
//                 reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//
//
//                 // memcpy(value.vec.data(), buf.data(),
//                 shape[0]*sizeof(type));
//
//                 // int result =
//                 detail::npy_api::get().PyArray_CopyInto_(ref.ptr(),
//                 buf.ptr());
//
//                 // if (result < 0)   // Copy failed!
//                 // {
//                    // PyErr_Clear();
//                    // return false;
//                 // }
//
//                 return true;
//             }
//
//             //Conversion part 2 (C++ -> Python)
//             static py::handle cast(const Vector<type>& src,
//             py::return_value_policy policy, py::handle parent)
//             {
//                 std::vector<size_t> shape(1);
//                 shape[0] = src.vec.size();
//
//                 std::vector<size_t> strides(1);
//                 strides[0] = sizeof(Real);
//
//                 py::array arr (std::move(shape), std::move(strides),
//                 src.dat);
//
//                 return arr.release();
//             }
//         };
//
//
//         template<typename type>
//         struct type_caster<Matrix<type>>
//         {
//             PYBIND11_TYPE_CASTER(Matrix<type>, _("Matrix<type>"));
//
//             // Conversion part 1 (Python -> C++)
//             bool load(py::handle src, bool convert)
//             {
//                 // if ( !convert and !py::array_t<T>::check_(src) )
//                     // return false;
//
//                 // auto buf = py::array_t<Real, py::array::c_style |
//                 py::array::forcecast>::ensure(src); auto buf =
//                 array::ensure(src);
//
//                 if ( !buf )
//                     return false;
//
//                 auto dims = buf.ndim();
//                 if ( dims != 2  )
//                     return false;
//
//                 std::vector<size_t> shape(2);
//                 shape[0] = buf.shape()[0];
//                 shape[1] = buf.shape()[1];
//
//                 value = Matrix<type> (shape[0], shape[1]);
//
//                 // auto ref =
//                 reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                 memcpy(value.vec.data(), buf.data(),
//                 shape[0]*shape[1]*sizeof(type));
//
//                 // int result =
//                 detail::npy_api::get().PyArray_CopyInto_(ref.ptr(),
//                 buf.ptr());
//
//                 // if (result < 0)   // Copy failed!
//                 // {
//                    // PyErr_Clear();
//                    // return false;
//                 // }
//
//                 return true;
//             }
//
//             //Conversion part 2 (C++ -> Python)
//             static py::handle cast(const Matrix<type>& src,
//             py::return_value_policy policy, py::handle parent)
//             {
//                 std::vector<size_t> shape(2);
//                 shape[0] = src.nrows;
//                 shape[1] = src.ncols;
//
//                 std::vector<size_t> strides(2);
//                 strides[0] = sizeof(Real)*shape[1];
//                 strides[1] = sizeof(Real);
//
//                 py::array arr (std::move(shape), std::move(strides),
//                 src.vec.data());
//
//                 return arr.release();
//             }
//         };
//     }
// }
