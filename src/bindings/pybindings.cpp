#include <iterator>

#include "../configure.hpp"   // ../ is required!
#include "tools/types.hpp"
#include "model/parameters/parameters.hpp"
#include "io/cpp/io_cpp_text.hpp"
#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/numpy.h"
#include "pybind11/eigen.h"
#include "pybind11/stl.h"
namespace py = pybind11;


PYBIND11_MAKE_OPAQUE (vector<LineProducingSpecies>);
PYBIND11_MAKE_OPAQUE (vector<CollisionPartner>);
// PYBIND11_MAKE_OPAQUE (vector<Matrix<Real>>);


PYBIND11_MODULE (core, module)
{
    // Module docstring
    module.doc() = "Core module of Magritte: a modern software library for 3D radiative transfer.";

    module.def(    "n_threads_avail", &paracabs::multi_threading::    n_threads_avail);
    module.def("set_n_threads_avail", &paracabs::multi_threading::set_n_threads_avail);

    // Define vector types
    py::bind_vector<vector<LineProducingSpecies>> (module, "vLineProducingSpecies");
    py::bind_vector<vector<CollisionPartner>>     (module, "vCollisionPartner");

    // Constants
    module.attr("CC")    = CC;
    module.attr("HH")    = HH;
    module.attr("KB")    = KB;
    module.attr("AMU")   = AMU;
    module.attr("T_CMB") = T_CMB;

    // Io, base class
    py::class_<Io> (module, "Io");


    // IoText
    py::class_<IoText, Io> (module, "IoText")
        // attributes
        .def_readonly ("io_file", &IoText::io_file)
        // constructor
        .def (py::init<const string &>());


    #if (PYTHON_IO)
        // IoPython
        py::class_<IoPython, Io> (module, "IoPython")
            // attributes
            .def_readonly ("implementation", &IoPython::implementation)
            .def_readonly ("io_file",        &IoPython::io_file)
            // constructor
            .def (py::init<const string &, const string &>());
    #endif

    // Solver
    py::class_<Solver> (module, "Solver")
        // attributes
        // functions
        // .def ("trace", &Solver::trace)
        // constructor
        // .def (py::init<const Size&, const Size&, const Size&>());
        .def (py::init<>());

    // Image
    py::class_<Image> (module, "Image")
        // attributes
        // functions
        .def_readonly  ("ray_nr", &Image::ray_nr)
        .def_readonly  ("ImX",    &Image::ImX)
        .def_readonly  ("ImY",    &Image::ImY)
        .def_readonly  ("I",      &Image::I)
        // constructor
        .def (py::init<const Geometry&, const Size&>());

    // Model
    py::class_<Model> (module, "Model")
        // attributes
        .def_readwrite ("parameters",     &Model::parameters)
        .def_readwrite ("geometry",       &Model::geometry)
        .def_readwrite ("chemistry",      &Model::chemistry)
        .def_readwrite ("lines",          &Model::lines)
        .def_readwrite ("a",          &Model::a)
        .def_readwrite ("b",          &Model::b)
        .def_readwrite ("c",          &Model::c)
        .def ("set",                  &Model::set)
        .def ("add",                  &Model::add)
        .def_readwrite ("thermodynamics", &Model::thermodynamics)
        .def_readwrite ("radiation",      &Model::radiation)
        // .def_readwrite ("neighbors_lists",&Model::neighbors_lists)
        .def_readonly  ("error_mean",     &Model::error_mean)
        .def_readonly  ("error_max",      &Model::error_max)
        .def_readonly  ("images",         &Model::images)
        // io (void (Pet::*)(int))
        .def ("read",  (void (Model::*)(void))            &Model::read )
        .def ("write", (void (Model::*)(void) const)      &Model::write)
        .def ("read",  (void (Model::*)(const Io&))       &Model::read )
        .def ("write", (void (Model::*)(const Io&) const) &Model::write)
        .def ("compute_inverse_line_widths",                                        &Model::compute_inverse_line_widths)
        .def ("compute_spectral_discretisation", (int (Model::*)(void            )) &Model::compute_spectral_discretisation)
        .def ("compute_spectral_discretisation", (int (Model::*)(const Real width)) &Model::compute_spectral_discretisation)
        .def ("compute_spectral_discretisation", (int (Model::*)(const long double nu_min, const long double nu_max)) &Model::compute_spectral_discretisation)
        .def ("compute_LTE_level_populations",                                      &Model::compute_LTE_level_populations)
        // .def ("compute_radiation_field",                                            &Model::compute_radiation_field)
        .def ("compute_radiation_field_feautrier_order_2",                          &Model::compute_radiation_field_feautrier_order_2)
        .def ("compute_radiation_field_shortchar_order_0",                          &Model::compute_radiation_field_shortchar_order_0)
        .def ("compute_Jeff",                                                       &Model::compute_Jeff)
        .def ("compute_level_populations_from_stateq",                              &Model::compute_level_populations_from_stateq)
        .def ("compute_level_populations",                                          &Model::compute_level_populations)
        .def ("compute_image",                                                      &Model::compute_image)
        .def ("set_eta_and_chi",                                                    &Model::set_eta_and_chi)
        .def ("set_boundary_condition",                                             &Model::set_boundary_condition)
        .def_readwrite ("eta",                &Model::eta)
        .def_readwrite ("chi",                &Model::chi)
        .def_readwrite ("boundary_condition", &Model::boundary_condition)
        // constructor
        .def (py::init<const string>())
        .def (py::init<>());

    // Parameters
    py::class_<Parameters> (module, "Parameters")
        // io
        .def_readwrite ("n_off_diag",         &Parameters::n_off_diag)
        .def_readwrite ("max_width_fraction", &Parameters::max_width_fraction)
        // setters
        .def ("set_model_name",               &Parameters::set_model_name          )
        .def ("set_dimension",                &Parameters::set_dimension           )
        .def ("set_npoints",                  &Parameters::set_npoints             )
        .def ("set_nrays",                    &Parameters::set_nrays               )
        .def ("set_hnrays",                   &Parameters::set_hnrays              )
        .def ("set_nrays_red",                &Parameters::set_nrays_red           )
        .def ("set_order_min",                &Parameters::set_order_min           )
        .def ("set_order_max",                &Parameters::set_order_max           )
        .def ("set_nboundary",                &Parameters::set_nboundary           )
        .def ("set_nfreqs",                   &Parameters::set_nfreqs              )
        .def ("set_nspecs",                   &Parameters::set_nspecs              )
        .def ("set_nlspecs",                  &Parameters::set_nlspecs             )
        .def ("set_nlines",                   &Parameters::set_nlines              )
        .def ("set_nquads",                   &Parameters::set_nquads              )
        .def ("set_pop_prec",                 &Parameters::set_pop_prec            )
        .def ("set_use_scattering",           &Parameters::set_use_scattering      )
        .def ("set_spherical_symmetry",       &Parameters::set_spherical_symmetry  )
        .def ("set_adaptive_ray_tracing",     &Parameters::set_adaptive_ray_tracing)
        // getters
        .def ("dimension",                    &Parameters::dimension           )
        .def ("npoints",                      &Parameters::npoints             )
        .def ("nrays",                        &Parameters::nrays               )
        .def ("hnrays",                       &Parameters::hnrays              )
        .def ("nrays_red",                    &Parameters::nrays_red           )
        .def ("order_min",                    &Parameters::order_min           )
        .def ("order_max",                    &Parameters::order_max           )
        .def ("nboundary",                    &Parameters::nboundary           )
        .def ("nfreqs",                       &Parameters::nfreqs              )
        .def ("nspecs",                       &Parameters::nspecs              )
        .def ("nlspecs",                      &Parameters::nlspecs             )
        .def ("nlines",                       &Parameters::nlines              )
        .def ("nquads",                       &Parameters::nquads              )
        .def ("pop_prec",                     &Parameters::pop_prec            )
        .def ("use_scattering",               &Parameters::use_scattering      )
        .def ("spherical_symmetry",           &Parameters::spherical_symmetry  )
        .def ("adaptive_ray_tracing",         &Parameters::adaptive_ray_tracing)
        // functions
        .def ("read",                         &Parameters::read                )
        .def ("write",                        &Parameters::write               )
        // constructor
        .def (py::init());


    // Geometry
    py::class_<Geometry> (module, "Geometry")
        // attributes
        .def_readwrite ("points",   &Geometry::points)
        .def_readwrite ("rays",     &Geometry::rays)
        .def_readwrite ("boundary", &Geometry::boundary)
        .def_readwrite ("lengths",  &Geometry::lengths)
        // io
        .def ("read",               &Geometry::read)
        .def ("write",              &Geometry::write)
        // functions
        // .def ("get_ray_lengths",     &Geometry::get_ray_lengths)
        // .def ("get_ray_lengths_gpu", &Geometry::get_ray_lengths_gpu)
        // constructor
        .def (py::init<>());


    // Points
    py::class_<Points> (module, "Points")
        // attributes
        .def_readwrite ("position",    &Points::position)
        .def_readwrite ("velocity",    &Points::velocity)
        .def_readwrite ("curr_neighbors", &Points::curr_neighbors)
        .def_readwrite ("multiscale", &Points::multiscale)
//@Frederik: I'm commenting these out for now :
// We should replace these with some methods from the Neighbors struct
//        .def_readwrite ("n_neighbors", &Points::n_neighbors)
//        .def_readwrite ("neighbors",   &Points::neighbors)
//        .def_readwrite ("nbs",         &Points::nbs)
        // io
        .def ("read",                  &Points::read)
        .def ("write",                 &Points::write)
        // constructor
        .def (py::init<>());


    // Neighbors
    py::class_<Neighbors> (module, "Neighbors")
        // attributes
        .def_readwrite ("n_neighbors",    &Neighbors::n_neighbors)
        .def_readwrite ("neighbors",    &Neighbors::neighbors)
        // functions
        .def("get_neighbors", &Neighbors::get_neighbors)
            // constructor
        .def (py::init<>());

    // Multiscale
    py::class_<Multiscale> (module, "Multiscale")
        // functions
        .def("set_all_neighbors", &Multiscale::set_all_neighbors)
        //.def("get_neighbors", &Multiscale::get_neighbors) TODO: figure out how to handle overloaded functions
        // constructor
        .def (py::init<>());


    // Rays
    py::class_<Rays> (module, "Rays")
        // attributes
        .def_readwrite ("direction", &Rays::direction)
        .def_readwrite ("antipod",   &Rays::antipod)
        .def_readwrite ("weight",    &Rays::weight)
        .def ("print",    &Rays::print)
        // io
        .def ("read",                &Rays::read)
        .def ("write",               &Rays::write)
        // constructor
        .def (py::init<>());


    // Boundary Condition
    py::enum_<BoundaryCondition>(module, "BoundaryCondition")
        .value("Zero",    Zero)
        .value("Thermal", Thermal)
        .value("CMB",     CMB)
        .export_values();


    // Boundary
    py::class_<Boundary> (module, "Boundary")
        // attributes
        .def_readwrite ("boundary2point",       &Boundary::boundary2point)
        .def_readwrite ("point2boundary",       &Boundary::point2boundary)
        .def_readwrite ("boundary_temperature", &Boundary::boundary_temperature)
        // functions
        .def ("set_boundary_condition",         &Boundary::set_boundary_condition)
        .def ("get_boundary_condition",         &Boundary::get_boundary_condition)
        // io
        .def ("read",                           &Boundary::read)
        .def ("write",                          &Boundary::write)
        // constructor
        .def (py::init<>());


    // Thermodynamics
    py::class_<Thermodynamics> (module, "Thermodynamics")
        // attributes
        .def_readwrite ("temperature", &Thermodynamics::temperature)
        .def_readwrite ("turbulence",  &Thermodynamics::turbulence)
        // io
        .def ("read",                  &Thermodynamics::read)
        .def ("write",                 &Thermodynamics::write)
        // constructor
        .def (py::init());


    // Temperature
    py::class_<Temperature> (module, "Temperature")
        // attributes
        .def_readwrite ("gas", &Temperature::gas)
        .def ("print",         &Temperature::print)
        // functions
        .def ("read",          &Temperature::read)
        .def ("write",         &Temperature::write)
        // constructor
        .def (py::init());


    // Turbulence
    py::class_<Turbulence> (module, "Turbulence")
        // attributes
        .def_readwrite ("vturb2", &Turbulence::vturb2)
        // functions
        .def ("read",             &Turbulence::read)
        .def ("write",            &Turbulence::write)
        // constructor
        .def (py::init());


    // Chemistry
    py::class_<Chemistry> (module, "Chemistry")
        // attributes
        .def_readwrite ("species", &Chemistry::species)
        // functions
        .def ("read",              &Chemistry::read)
        .def ("write",             &Chemistry::write)
        // constructor
        .def (py::init());


    // Species
    py::class_<Species> (module, "Species")
        // attributes
        .def_readwrite ("symbol",    &Species::symbol)
        .def_readwrite ("abundance", &Species::abundance)
        // functions
        .def ("read",                &Species::read)
        .def ("write",               &Species::write)
        // constructor
        .def (py::init());


    // Lines
    py::class_<Lines> (module, "Lines")
        // attributes
        .def_readwrite ("lineProducingSpecies", &Lines::lineProducingSpecies)
        .def_readwrite ("emissivity",           &Lines::emissivity)
        .def_readwrite ("opacity",              &Lines::opacity)
        .def_readwrite ("inverse_width",        &Lines::inverse_width)
        .def_readwrite ("line",                 &Lines::line)
        // functions
        .def ("read",                           &Lines::read)
        .def ("write",                          &Lines::write)
        .def ("set_emissivity_and_opacity",     &Lines::set_emissivity_and_opacity)
        // constructor
        .def (py::init<>());


    // LineProducingSpecies
    py::class_<LineProducingSpecies> (module, "LineProducingSpecies")
        // attributes
        .def_readwrite ("linedata",         &LineProducingSpecies::linedata)
        .def_readwrite ("quadrature",       &LineProducingSpecies::quadrature)
        .def_readwrite ("Lambda",           &LineProducingSpecies::lambda) // "lambda" is invalid in Python, use "Lambda"
        .def_readwrite ("Jeff",             &LineProducingSpecies::Jeff)
        .def_readwrite ("Jdif",             &LineProducingSpecies::Jdif)
        .def_readwrite ("Jlin",             &LineProducingSpecies::Jlin)
        .def_readwrite ("nr_line",          &LineProducingSpecies::nr_line)
        .def_readwrite ("population",       &LineProducingSpecies::population)
        .def_readwrite ("population_tot",   &LineProducingSpecies::population_tot)
        .def_readwrite ("population_prev1", &LineProducingSpecies::population_prev1)
        .def_readwrite ("population_prev2", &LineProducingSpecies::population_prev2)
        .def_readwrite ("population_prev3", &LineProducingSpecies::population_prev3)
        .def_readwrite ("populations",      &LineProducingSpecies::populations)
        .def_readwrite ("RT",               &LineProducingSpecies::RT)
        .def_readwrite ("LambdaStar",       &LineProducingSpecies::LambdaStar)
        .def_readwrite ("LambdaTest",       &LineProducingSpecies::LambdaTest)
        // functions
        .def ("read",                       &LineProducingSpecies::read)
        .def ("write",                      &LineProducingSpecies::write)
        .def ("index",                      &LineProducingSpecies::index)
        // constructor
        .def (py::init<>());


    // Lambda
    py::class_<Lambda> (module, "Lambda")
        // attributes
        .def_readwrite ("Ls",   &Lambda::Ls)
        .def_readwrite ("nr",   &Lambda::nr)
        .def_readwrite ("size", &Lambda::size)
        .def_readwrite ("Lss",  &Lambda::Lss)
        .def_readwrite ("nrs",  &Lambda::nrs)
        // functions
        .def ("add_element",    &Lambda::add_element)
        .def ("linearize_data", &Lambda::linearize_data)
        .def ("MPI_gather",     &Lambda::MPI_gather)
        // constructor
        .def (py::init<>());


    // Quadrature
    py::class_<Quadrature> (module, "Quadrature")
        // attributes
        .def_readwrite ("roots",   &Quadrature::roots)
        .def_readwrite ("weights", &Quadrature::weights)
        // functions
        .def ("read",              &Quadrature::read)
        .def ("write",             &Quadrature::write)
        // constructor
        .def (py::init<>());


    // Linedata
    py::class_<Linedata> (module, "Linedata")
        // attributes
        .def_readwrite ("num",          &Linedata::num)
        .def_readwrite ("sym",          &Linedata::sym)
        .def_readwrite ("inverse_mass", &Linedata::inverse_mass)
        .def_readwrite ("nlev",         &Linedata::nlev)
        .def_readwrite ("nrad",         &Linedata::nrad)
        .def_readwrite ("irad",         &Linedata::irad)
        .def_readwrite ("jrad",         &Linedata::jrad)
        .def_readwrite ("energy",       &Linedata::energy)
        .def_readwrite ("weight",       &Linedata::weight)
        .def_readwrite ("frequency",    &Linedata::frequency)
        .def_readwrite ("A",            &Linedata::A)
        .def_readwrite ("Ba",           &Linedata::Ba)
        .def_readwrite ("Bs",           &Linedata::Bs)
        .def_readwrite ("ncolpar",      &Linedata::ncolpar)
        .def_readwrite ("colpar",       &Linedata::colpar)
        // functions
        .def ("read",                   &Linedata::read)
        .def ("write",                  &Linedata::write)
        // constructor
        .def (py::init<>());


    // Colpartner
    py::class_<CollisionPartner> (module, "CollisionPartner")
        // attributes
        .def_readwrite ("num_col_partner", &CollisionPartner::num_col_partner)
        .def_readwrite ("orth_or_para_H2", &CollisionPartner::orth_or_para_H2)
        .def_readwrite ("ntmp",            &CollisionPartner::ntmp)
        .def_readwrite ("ncol",            &CollisionPartner::ncol)
        .def_readwrite ("icol",            &CollisionPartner::icol)
        .def_readwrite ("jcol",            &CollisionPartner::jcol)
        .def_readwrite ("tmp",             &CollisionPartner::tmp)
        .def_readwrite ("Ce",              &CollisionPartner::Ce)
        .def_readwrite ("Cd",              &CollisionPartner::Cd)
        // functions
        .def ("read",                      &CollisionPartner::read)
        .def ("write",                     &CollisionPartner::write)
        // constructor
        .def (py::init<>());


    // Radiation
    py::class_<Radiation> (module, "Radiation")
        // attributes
        .def_readwrite ("frequencies", &Radiation::frequencies)
        // .def_readwrite ("u",           &Radiation::u)
        // .def_readwrite ("v",           &Radiation::v)
        .def_readwrite ("I_bdy",       &Radiation::I_bdy)
        .def_readwrite ("I",           &Radiation::I)
        .def_readwrite ("u",           &Radiation::u)
        .def_readwrite ("v",           &Radiation::v)
        .def_readwrite ("J",           &Radiation::J)
        // functions
        .def ("read",                  &Radiation::read)
        .def ("write",                 &Radiation::write)
        // constructor
        .def (py::init());


    // Frequencies
    py::class_<Frequencies> (module, "Frequencies")
        // attributes
        .def_readwrite ("nu", &Frequencies::nu)
        // functions
        .def ("read",         &Frequencies::read)
        .def ("write",        &Frequencies::write)
        // constructor
        .def (py::init());


    // Vector <Size>
    py::class_<Vector<Size>> (module, "VSize", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Vector<Size> &v) -> py::buffer_info
            {
                return py::buffer_info(
                    v.vec.data(),                          // Pointer to buffer
                    sizeof(Size),                          // Size of one element
                    py::format_descriptor<Size>::format(), // Python struct-style format descriptor
                    1,                                     // Number of dimensions
                    {v.vec.size()},                        // Buffer dimensions
                    {sizeof(Size)}                         // Strides (in bytes) for each index
                );
            }
        )
        // functions
        .def ("set", &Vector<Size>::set_1D_array)
        // constructor
        .def (py::init());


    // Vector <Real>
    py::class_<Vector<Real>> (module, "VReal", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Vector<Real> &v) -> py::buffer_info
            {
                return py::buffer_info(
                    v.vec.data(),                          // Pointer to buffer
                    sizeof(Real),                          // Size of one element
                    py::format_descriptor<Real>::format(), // Python struct-style format descriptor
                    1,                                     // Number of dimensions
                    {v.vec.size()},                        // Buffer dimensions
                    {sizeof(Real)}                         // Strides (in bytes) for each index
                );
            }
        )
        // functions
        .def ("set", &Vector<Real>::set_1D_array)
        // constructor
        .def (py::init());


    // Matrix <Real>
    py::class_<Matrix<Real>, Vector<Real>> (module, "MReal", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Matrix<Real> &m) -> py::buffer_info
            {
                return py::buffer_info(
                    m.vec.data(),                                                 // Pointer to buffer
                    sizeof(Real),                                                 // Size of one element
                    py::format_descriptor<Real>::format(),                        // Python struct-style format descriptor
                    2,                                                            // Number of dimensions
                    py::detail::any_container<ssize_t>({m.nrows,
                                                        m.ncols}),                // Buffer dimensions
                    py::detail::any_container<ssize_t>({sizeof(Real)*m.ncols,
                                                        sizeof(Real)         })   // Strides (in bytes) for each index
                );
            }
        )
        .def_readwrite ("vec",   &Vector<Real>::vec)
        .def_readwrite ("nrows", &Matrix<Real>::nrows)
        .def_readwrite ("ncols", &Matrix<Real>::ncols)
        // functions
        .def ("set", &Matrix<Real>::set_2D_array)
        // constructor
        .def (py::init());


    // Tensor <Real>
    py::class_<Tensor<Real>, Vector<Real>> (module, "TReal", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Tensor<Real> &t) -> py::buffer_info
            {
                return py::buffer_info(
                    t.vec.data(),                                                        // Pointer to buffer
                    sizeof(Real),                                                        // Size of one element
                    py::format_descriptor<Real>::format(),                               // Python struct-style format descriptor
                    3,                                                                   // Number of dimensions
                    py::detail::any_container<ssize_t>({t.nrows,
                                                        t.ncols,
                                                        t.depth }),                      // Buffer dimensions
                    py::detail::any_container<ssize_t>({sizeof(Real)*t.ncols*t.depth,
                                                        sizeof(Real)*t.depth,
                                                        sizeof(Real)                 })   // Strides (in bytes) for each index
                );
            }
        )
        .def_readwrite ("vec",   &Vector<Real>::vec)
        .def_readwrite ("nrows", &Tensor<Real>::nrows)
        .def_readwrite ("ncols", &Tensor<Real>::ncols)
        .def_readwrite ("depth", &Tensor<Real>::depth)
        // functions
        .def ("set", &Tensor<Real>::set_3D_array)
        // constructor
        .def (py::init());


    // Vector <Vector3D>
    py::class_<Vector<Vector3D>> (module, "VVector3D", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Vector<Vector3D> &v) -> py::buffer_info
            {
                return py::buffer_info(
                    v.vec.data(),                                             // Pointer to buffer
                    sizeof(double),                                           // Size of one element
                    py::format_descriptor<double>::format(),                  // Python struct-style format descriptor
                    2,                                                        // Number of dimensions
                    py::detail::any_container<ssize_t>({v.vec.size(),
                                                        3            }),      // Buffer dimensions
                    py::detail::any_container<ssize_t>({sizeof(double)*3,
                                                        sizeof(double)   })   // Strides (in bytes) for each index
                );
            }
        )
        // functions
        .def ("set", (void (Vector<Vector3D>::*)(py::array_t<double, py::array::c_style | py::array::forcecast>)) &Vector<Vector3D>::set_2D_array)
        // constructor
        .def (py::init());
}


//namespace pybind11
//{
//    namespace detail
//    {
//        template<> struct type_caster<Vector3D>
//        {
//            PYBIND11_TYPE_CASTER(Vector3D, _("Vector3D"));
//
//            // Conversion part 1 (Python -> C++)
//            bool load(py::handle src, bool convert)
//            {
//                // if ( !convert and !py::array_t<T>::check_(src) )
//                    // return false;
//
//                // auto buf = py::array_t<Real, py::array::c_style | py::array::forcecast>::ensure(src);
//                auto buf = array::ensure(src);
//
//                if ( !buf )
//                    return false;
//
//                auto dims = buf.ndim();
//                if ( dims != 1  )
//                    return false;
//
//                // std::vector<size_t> shape(1);
//                // shape[0] = buf.shape()[0];
//                Real* data = (Real*) buf.data();
//
//                value = Vector3D (data[0], data[1], data[2]);
//
//                // auto ref = reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                // memcpy(value.vec.data(), buf.data(), shape[0]*sizeof(Vector3D));
//
//                // int result = detail::npy_api::get().PyArray_CopyInto_(ref.ptr(), buf.ptr());
//
//                // if (result < 0)   // Copy failed!
//                // {
//                   // PyErr_Clear();
//                   // return false;
//                // }
//
//                return true;
//            }
//
//            //Conversion part 2 (C++ -> Python)
//            static py::handle cast(const Vector3D& src, py::return_value_policy policy, py::handle parent)
//            {
//                std::vector<size_t> shape(1);
//                shape[0] = 3;
//
//                std::vector<size_t> strides(1);
//                strides[0] = sizeof(Real);
//
//                py::array arr (std::move(shape), std::move(strides), &(src.data[0]));
//
//                return arr.release();
//            }
//        };
//
//
//        template<> struct type_caster<Vector<Vector3D>>
//        {
//            PYBIND11_TYPE_CASTER(Vector<Vector3D>, _("Vector<Vector3D>"));
//
//            // Conversion part 1 (Python -> C++)
//            bool load(py::handle src, bool convert)
//            {
//                auto buf = array::ensure(src);
//                if (!buf) return false;
//
//                auto dims = buf.ndim();
//                if (dims != 2) return false;
//
//                std::vector<size_t> shape(2);
//                shape[0] = buf.shape()[0];
//                shape[1] = buf.shape()[1];
//
//                if (shape[1] != 3) return false;
//
//                // value = Vector<Vector3D> (shape[0]);
//                value.resize(shape[0]);
//
//                // Vector3D>* data = (Vector<Vector3D>*) buf.data();
//
//                // for (size_t i = 0; i < shape[0]; i++)
//                // {
//                    // value[i] = Vector3D()
//                // }
//
//                // Vector<Vector3D>* data = (Vector<Vector3D>*) buf.data();
//
//                // auto ref = py::reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                memcpy(&(value.vec[0].data[0]), buf.data(), shape[0]*shape[1]*sizeof(Real));
//
//                // int result = detail::npy_api::get().PyArray_CopyInto_(ref.ptr(), buf.ptr());
//
//                // if (result < 0)   // Copy failed!
//                // {
//                   // PyErr_Clear();
//                   // return false;
//                // }
//
//                return true;
//            }
//
//            //Conversion part 2 (C++ -> Python)
//            static py::handle cast(const Vector<Vector3D>& src, py::return_value_policy policy, py::handle parent)
//            {
//                std::vector<size_t> shape(2);
//                shape[0] = src.vec.size();
//                shape[1] = 3;
//
//                std::vector<size_t> strides(2);
//                strides[0] = sizeof(Real)*shape[1];
//                strides[1] = sizeof(Real);
//
//                py::array arr (std::move(shape), std::move(strides), &(src.vec[0].data[0]));
//
//                return arr.release();
//            }
//        };
//
//
//        template<typename type>
//        struct type_caster<Vector<type>>
//        {
//            PYBIND11_TYPE_CASTER(Vector<type>, _("Vector<type>"));
//
//            // Conversion part 1 (Python -> C++)
//            bool load(py::handle src, bool convert)
//            {
//                // if ( !convert and !py::array_t<T>::check_(src) )
//                    // return false;
//
//                // auto buf = py::array_t<Real, py::array::c_style | py::array::forcecast>::ensure(src);
//                auto buf = py::array_t<type, py::array::c_style | py::array::forcecast>::ensure(src);
//
//                if ( !buf )
//                    return false;
//
//                auto dims = buf.ndim();
//                if ( dims != 1  )
//                    return false;
//
//                std::vector<size_t> shape(1);
//                shape[0] = buf.shape()[0];
//
//                cout << "----------- value.data() IN  " << value.vec.data() << endl;
//                cout << "----------- value.dat    IN  " << value.dat        << endl;
//
//                value = Vector<type> (shape[0]);
//
//                // std::copy (buf.data(), buf.data()+buf.size(), value.vec.data());
//                cout << "buf.size() = " << buf.size() << endl;
//                cout << "shape[0]   = " << shape[0]   << endl;
//                cout << "----------- value.data() IN  " << value.vec.data() << endl;
//                cout << "----------- value.dat    IN  " << value.dat        << endl;
//                value.set_dat();
//                cout << "----------- value.data() IN  " << value.vec.data() << endl;
//                cout << "----------- value.dat    IN  " << value.dat        << endl;
//                std::copy (buf.data(), buf.data()+buf.size(), value.dat);
//                cout << "----------- value.data() OUT " << value.vec.data() << endl;
//                cout << "----------- value.dat    OUT " << value.dat        << endl;
//
//                // auto ref = reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//
//
//                // memcpy(value.vec.data(), buf.data(), shape[0]*sizeof(type));
//
//                // int result = detail::npy_api::get().PyArray_CopyInto_(ref.ptr(), buf.ptr());
//
//                // if (result < 0)   // Copy failed!
//                // {
//                   // PyErr_Clear();
//                   // return false;
//                // }
//
//                return true;
//            }
//
//            //Conversion part 2 (C++ -> Python)
//            static py::handle cast(const Vector<type>& src, py::return_value_policy policy, py::handle parent)
//            {
//                std::vector<size_t> shape(1);
//                shape[0] = src.vec.size();
//
//                std::vector<size_t> strides(1);
//                strides[0] = sizeof(Real);
//
//                py::array arr (std::move(shape), std::move(strides), src.dat);
//
//                return arr.release();
//            }
//        };
//
//
//        template<typename type>
//        struct type_caster<Matrix<type>>
//        {
//            PYBIND11_TYPE_CASTER(Matrix<type>, _("Matrix<type>"));
//
//            // Conversion part 1 (Python -> C++)
//            bool load(py::handle src, bool convert)
//            {
//                // if ( !convert and !py::array_t<T>::check_(src) )
//                    // return false;
//
//                // auto buf = py::array_t<Real, py::array::c_style | py::array::forcecast>::ensure(src);
//                auto buf = array::ensure(src);
//
//                if ( !buf )
//                    return false;
//
//                auto dims = buf.ndim();
//                if ( dims != 2  )
//                    return false;
//
//                std::vector<size_t> shape(2);
//                shape[0] = buf.shape()[0];
//                shape[1] = buf.shape()[1];
//
//                value = Matrix<type> (shape[0], shape[1]);
//
//                // auto ref = reinterpret_steal<array>(eigen_ref_array<props>(value));
//
//                memcpy(value.vec.data(), buf.data(), shape[0]*shape[1]*sizeof(type));
//
//                // int result = detail::npy_api::get().PyArray_CopyInto_(ref.ptr(), buf.ptr());
//
//                // if (result < 0)   // Copy failed!
//                // {
//                   // PyErr_Clear();
//                   // return false;
//                // }
//
//                return true;
//            }
//
//            //Conversion part 2 (C++ -> Python)
//            static py::handle cast(const Matrix<type>& src, py::return_value_policy policy, py::handle parent)
//            {
//                std::vector<size_t> shape(2);
//                shape[0] = src.nrows;
//                shape[1] = src.ncols;
//
//                std::vector<size_t> strides(2);
//                strides[0] = sizeof(Real)*shape[1];
//                strides[1] = sizeof(Real);
//
//                py::array arr (std::move(shape), std::move(strides), src.vec.data());
//
//                return arr.release();
//            }
//        };
//    }
//}
