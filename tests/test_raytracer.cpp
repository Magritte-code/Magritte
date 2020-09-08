#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"


//__global__
//void kernel ()
//{
//    return;
//}

//inline void test (Geometry &geometry)
//{
//    geometry.points.position.vec[0].print();
//    geometry.points.position.vec[1].print();
//
//    Geometry* geometry_copy = (Geometry*) pc::accelerator::malloc(sizeof(Geometry));
//    pc::accelerator::memcpy_to_accelerator (geometry_copy, &geometry, sizeof(Geometry));
//    kernel<<<1,1>>> (geometry_copy);
//
//    points.position.vec[0].print();
//    points.position.copy_ptr_to_vec ();
//    points.position.vec[0].print();
//}


int main ()
{
    cout << "Running test_raytracing..." << endl;

    pc::accelerator::list_accelerators();

    IoPython io ("hdf5", "/home/frederik/Magritte_all/Models/Benchmarks/5_GPU_ray_tracer/model.hdf5");


    cout << "sizeof Model    = " << sizeof(Model)    << endl;
    cout << "sizeof geometry = " << sizeof(Geometry) << endl;
    cout << "sizeof points   = " << sizeof(Points)   << endl;
    cout << "sizeof Vector3D = " << sizeof(Vector3D) << endl;

//    model.read(io);

//    Vector<double> t;

//    t.resize (100);

    cout << "n threads = " << paracabs::multi_threading::n_threads() << endl;

    Solver solver (10000, 100);

    Model model;
    model.read (io);

//    Parameters params1;
//    params1.a() = 5;

//    model.parameters.set_a(5);

    solver.trace (model);
//
//    for (Size i = 0; i < 10; i++)
//    {
//        cout << model.geometry.lengths[i] << endl;
//    }

//    Size1 lengths = model.geometry.get_ray_lengths ();
//    Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);

    for (Size i = 0; i < 100; i++)
    {
        cout << model.geometry.lengths[i] << endl;
    }

    cout << "Done." << endl;

//    cout << params1.a() << endl;

//    params1.a() = 8;

//    cout << params1.a() << endl;

//    Parameters params2;
//    cout << params2.a() << endl;

//    params2.a() = 13;

//    cout << params1.a() << endl;


    return (0);
}