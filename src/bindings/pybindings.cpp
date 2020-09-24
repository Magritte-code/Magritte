#include "../configure.hpp"   // ../ is required!
#include "tools/types.hpp"
#include "model/parameters/parameters.hpp"
#include "io/cpp/io_cpp_text.hpp"
#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/eigen.h"
#include "pybind11/stl.h"
namespace py = pybind11;


PYBIND11_MODULE (core, module)
{
    // Module docstring
    module.doc() = "Core module of Magritte: a modern software library for 3D radiative transfer.";


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
        .def ("trace", &Solver::trace)
        // constructor
        .def (py::init<const Size&, const Size&>());

    // Model
    py::class_<Model> (module, "Model")
        // attributes
        .def_readwrite ("parameters", &Model::parameters)
        .def_readwrite ("geometry",   &Model::geometry)
        // io
        .def ("read",                 &Model::read)
        .def ("write",                &Model::write)
        // constructor
        .def (py::init());

    // Parameters
    py::class_<Parameters> (module, "Parameters")
        // io
        .def_readwrite ("n_off_diag",         &Parameters::n_off_diag)
        .def_readwrite ("max_width_fraction", &Parameters::max_width_fraction)
        // setters
        .def ("set_npoints",                  &Parameters::set_npoints             )
        .def ("set_nrays",                    &Parameters::set_nrays               )
        .def ("set_nrays_red",                &Parameters::set_nrays_red           )
        .def ("set_order_min",                &Parameters::set_order_min           )
        .def ("set_order_max",                &Parameters::set_order_max           )
        .def ("set_nboundary",                &Parameters::set_nboundary           )
        .def ("set_nfreqs",                   &Parameters::set_nfreqs              )
        .def ("set_nfreqs_red",               &Parameters::set_nfreqs_red          )
        .def ("set_nspecs",                   &Parameters::set_nspecs              )
        .def ("set_nlspecs",                  &Parameters::set_nlspecs             )
        .def ("set_nlines",                   &Parameters::set_nlines              )
        .def ("set_nquads",                   &Parameters::set_nquads              )
        .def ("set_pop_prec",                 &Parameters::set_pop_prec            )
        .def ("set_use_scattering",           &Parameters::set_use_scattering      )
        .def ("set_spherical_symmetry",       &Parameters::set_spherical_symmetry  )
        .def ("set_adaptive_ray_tracing",     &Parameters::set_adaptive_ray_tracing)
        // getters
        .def ("npoints",                      &Parameters::npoints             )
        .def ("nrays",                        &Parameters::nrays               )
        .def ("nrays_red",                    &Parameters::nrays_red           )
        .def ("order_min",                    &Parameters::order_min           )
        .def ("order_max",                    &Parameters::order_max           )
        .def ("nboundary",                    &Parameters::nboundary           )
        .def ("nfreqs",                       &Parameters::nfreqs              )
        .def ("nfreqs_red",                   &Parameters::nfreqs_red          )
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
        .def ("get_ray_lengths",     &Geometry::get_ray_lengths)
        .def ("get_ray_lengths_gpu", &Geometry::get_ray_lengths_gpu)
        // constructor
        .def (py::init<Parameters&>());


    // Points
    py::class_<Points> (module, "Points")
        // attributes
        .def_readwrite ("position",    &Points::position)
        .def_readwrite ("velocity",    &Points::velocity)
//@Frederik: I'm commenting these out for now :
// We should replace these with some methods from the Neighbors struct
//        .def_readwrite ("n_neighbors", &Points::n_neighbors)
//        .def_readwrite ("neighbors",   &Points::neighbors)
//        .def_readwrite ("nbs",         &Points::nbs)
        // io
        .def ("read",                  &Points::read)
        .def ("write",                 &Points::write)
        // constructor
        .def (py::init<Parameters&>());


    // Rays
    py::class_<Rays> (module, "Rays")
        // attributes
        .def_readwrite ("direction", &Rays::direction)
        .def_readwrite ("antipod",   &Rays::antipod)
        .def_readwrite ("weight",    &Rays::weight)
        // io
        .def ("read",                &Rays::read)
        .def ("write",               &Rays::write)
        // constructor
        .def (py::init<Parameters&>());


    // Boundary
    py::class_<Boundary> (module, "Boundary")
        // attributes
        .def_readwrite ("boundary2point", &Boundary::boundary2point)
        .def_readwrite ("point2boundary", &Boundary::point2boundary)
//        .def_readwrite ("is_on_boundary", &Boundary::is_on_boundary)
        // io
        .def ("read",                     &Boundary::read)
        .def ("write",                    &Boundary::write)
        // constructor
        .def (py::init<Parameters&>());


    // Vector <Size>
    py::class_<Vector<Size>> (module, "Vector_Size", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Vector<Size> &v) -> py::buffer_info
            {
                return py::buffer_info(
                    v.dat,                                 // Pointer to buffer
                    sizeof(Size),                          // Size of one element
                    py::format_descriptor<Size>::format(), // Python struct-style format descriptor
                    1,                                     // Number of dimensions
                    {v.vec.size()},                        // Buffer dimensions
                    {sizeof(Size)}                         // Strides (in bytes) for each index
                );
            }
        )
        // functions
        .def ("resize",        &Vector<Size>::resize)
        // attributes
        .def_readwrite ("vec", &Vector<Size>::vec)
        // constructor
        .def (py::init());


    // Vector <Real>
    py::class_<Vector<Real>> (module, "Vector_Real", py::buffer_protocol())
        // buffer
        .def_buffer(
            [](Vector<Real> &v) -> py::buffer_info
            {
                return py::buffer_info(
                    v.dat,                                 // Pointer to buffer
                    sizeof(Real),                          // Size of one element
                    py::format_descriptor<Real>::format(), // Python struct-style format descriptor
                    1,                                     // Number of dimensions
                    {v.vec.size()},                        // Buffer dimensions
                    {sizeof(Real)}                         // Strides (in bytes) for each index
                );
            }
        )
        // functions
        .def ("resize",        &Vector<Real>::resize)
        // attributes
        .def_readwrite ("vec", &Vector<Real>::vec)
        // constructor
        .def (py::init());
}
