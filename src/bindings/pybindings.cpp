#include "../configure.hpp"   // ../ is required!
#include "tools/types.hpp"
#include "io/cpp/io_cpp_text.hpp"
#include "io/python/io_python.hpp"
#include "model/model.hpp"

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

    bool truth = PYTHON_IO;

        #if (PYTHON_IO)
        // IoPython
        py::class_<IoPython, Io> (module, "IoPython")
            // attributes
            .def_readonly ("implementation", &IoPython::implementation)
            .def_readonly ("io_file",        &IoPython::io_file)
            // constructor
            .def (py::init<const string &, const string &>());
        #endif


    // Model
    py::class_<Model> (module, "Model")
        // attributes
        .def_readwrite ("geometry", &Model::geometry)
        // constructor
        .def (py::init())
        // io
        .def ("read",               &Model::read)
        .def ("write",              &Model::write);


    // Geometry
    py::class_<Geometry> (module, "Geometry")
        // attributes
        .def_readwrite ("points",   &Geometry::points)
        .def_readwrite ("rays",     &Geometry::rays)
        .def_readwrite ("boundary", &Geometry::boundary)
        // io
        .def ("read",               &Geometry::read)
        .def ("write",              &Geometry::write)
        // functions
        .def ("get_ray_lengths",    &Geometry::get_ray_lengths)
        // constructor
        .def (py::init());


// Points
    py::class_<Points> (module, "Points")
        // attributes
        .def_readwrite ("position",    &Points::position)
        .def_readwrite ("velocity",    &Points::velocity)
        .def_readwrite ("n_neighbors", &Points::n_neighbors)
        .def_readwrite ("neighbors",   &Points::neighbors)
        // constructor
        .def (py::init())
        // functions
        .def ("read",                  &Points::read)
        .def ("write",                 &Points::write);


    // Rays
    py::class_<Rays> (module, "Rays")
        // attributes
        .def_readwrite ("direction", &Rays::direction)
        .def_readwrite ("antipod",   &Rays::antipod)
        .def_readwrite ("weight",    &Rays::weight)
        // constructor
        .def (py::init())
        // io
        .def ("read",                &Rays::read)
        .def ("write",               &Rays::write);


    // Boundary
    py::class_<Boundary> (module, "Boundary")
        // attributes
        .def_readwrite ("boundary2point", &Boundary::boundary2point)
        .def_readwrite ("point2boundary", &Boundary::point2boundary)
        .def_readwrite ("is_on_boundary", &Boundary::is_on_boundary)
        // constructor
        .def (py::init())
        // functions
        .def ("read",                     &Boundary::read)
        .def ("write",                    &Boundary::write);
}
