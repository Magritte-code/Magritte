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
        .def ("solve", &Solver::solve)
        // constructor
        .def (py::init<const Size&, const Size&>());

    // Model
    py::class_<Model> (module, "Model")
        // attributes
        .def_readwrite ("parameters",     &Model::parameters)
        .def_readwrite ("geometry",       &Model::geometry)
        .def_readwrite ("chemistry",      &Model::chemistry)
        .def_readwrite ("lines",          &Model::lines)
        .def_readwrite ("thermodynamics", &Model::thermodynamics)
        .def_readwrite ("radiation",      &Model::radiation)
        // io
        .def ("read",                     &Model::read)
        .def ("write",                    &Model::write)
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
        .def_readwrite ("n_neighbors", &Points::n_neighbors)
        .def_readwrite ("neighbors",   &Points::neighbors)
        .def_readwrite ("nbs",         &Points::nbs)
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
        .def_readwrite ("boundary_condition",   &Boundary::boundary_condition)
        .def_readwrite ("boundary_temperature", &Boundary::boundary_temperature)
        // io
        .def ("read",                           &Boundary::read)
        .def ("write",                          &Boundary::write)
        // constructor
        .def (py::init<Parameters&>());


    // Thermodynamics
    py::class_<Thermodynamics> (module, "Thermodynamics")
        // attributes
        .def_readwrite ("temperature", &Thermodynamics::temperature)
        .def_readwrite ("turbulence",  &Thermodynamics::turbulence)
        // functions
        .def ("read",                  &Thermodynamics::read)
        .def ("write",                 &Thermodynamics::write)
        // constructor
        .def (py::init());


    // Temperature
    py::class_<Temperature> (module, "Temperature")
        // attributes
        .def_readwrite ("gas", &Temperature::gas)
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
        .def_readwrite ("Jlin",             &LineProducingSpecies::Jlin)
        .def_readwrite ("nr_line",          &LineProducingSpecies::nr_line)
        .def_readwrite ("population",       &LineProducingSpecies::population)
        .def_readwrite ("population_tot",   &LineProducingSpecies::population_tot)
        .def_readwrite ("population_prev1", &LineProducingSpecies::population_prev1)
        .def_readwrite ("population_prev2", &LineProducingSpecies::population_prev2)
        .def_readwrite ("population_prev3", &LineProducingSpecies::population_prev3)
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
