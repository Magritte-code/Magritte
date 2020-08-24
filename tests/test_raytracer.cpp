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

    IoPython io ("hdf5", "/home/frederik/Magritte_all/Models/Benchmarks/5_GPU_ray_tracer/test_model_0.hdf5");


    cout << "sizeof Model    = " << sizeof(Model)    << endl;
    cout << "sizeof geometry = " << sizeof(Geometry) << endl;
    cout << "sizeof points   = " << sizeof(Points)   << endl;
    cout << "sizeof Vector3D = " << sizeof(Vector3D) << endl;

//    model.read(io);

//    Vector<double> t;

//    t.resize (100);


    Model model;
    model.read (io);
//    model.geometry.test();

//    Vector<Vector<size_t>> test (10, Vector<size_t> (20, 3));
//
//    test.copy_vec_to_ptr ();
//    for (auto& t : test.vec) {t.copy_vec_to_ptr ();}

//    Long2 lengths = model.geometry.get_ray_lengths ();
    Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);

    cout << "Done." << endl;

    return (0);
}