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


int main ()
{
    cout << "Running test_raytracing..." << endl;

    pc::accelerator::list_accelerators();

//    kernel<<<1,1>>>();
    IoPython io ("hdf5", "/home/frederik/Magritte_all/Models/Benchmarks/5_GPU_ray_tracer/test_model.hdf5");

    Model model;
    model.read(io);

    Long2 lengths = model.geometry.get_ray_lengths ();

    cout << "Done." << endl;

    return (0);
}