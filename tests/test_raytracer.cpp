#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"


int main (int argc, char **argv)
{
    cout << "Running test_raytracing..." << endl;

    /// Store model name
    const string modelName = argv[1];

    cout << modelName << endl;

    pc::accelerator::list_accelerators();

    IoPython io ("hdf5", modelName);


    cout << "sizeof Model    = " << sizeof(Model)    << endl;
    cout << "sizeof geometry = " << sizeof(Geometry) << endl;
    cout << "sizeof points   = " << sizeof(Points)   << endl;
    cout << "sizeof Vector3D = " << sizeof(Vector3D) << endl;

    cout << "n threads = " << paracabs::multi_threading::n_threads_avail() << endl;


    Model model;
    model.read (io);

    Solver solver (10000, 10000);
    solver.trace  (model);

    // Size1 lengths = model.geometry.get_ray_lengths ();
    // Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);

    for (Size i = 0; i < 100; i++)
    {
       cout << model.geometry.lengths[i] << endl;
        // cout << lengths[i] << endl;
    }

    cout << "Done." << endl;

    return (0);
}
