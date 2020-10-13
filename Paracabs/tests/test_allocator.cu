#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
using namespace paracabs::accelerator;
//using paracabs::datatypes::vector3d;
//using paracabs::datatypes::array1d;
using paracabs::datatypes::MemTypeDefault;
using paracabs::datatypes::MemTypeAccelerator;
using paracabs::datatypes::my_copy;

//template <typename T>
//typedef  vector<T, > Vector<T>;

//typedef vector3d<double, MemTypeDefault>     Vector3d;
//typedef vector3d<double, MemTypeAccelerator> Vector3d_accel;

//typedef array1d <Vector3d,       MemTypeDefault>     arrvec3d;
//typedef array1d <Vector3d_accel, MemTypeAccelerator> arrvec3d_accel;

//typedef array1d <double, MemTypeDefault>     arrdouble;
//typedef array1d <double, MemTypeAccelerator> arrdouble_accel;


__global__ void printKernel (vector <double, paracabs::allocator<double, MemTypeAccelerator>> vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        printf("test\n");
    }
}


int main ()
{
    cout << "Paracabs test allocator."     << endl;
    cout << "Number of GPUs = " << nGPUs() << endl;

    list_accelerators();

    const size_t size = 10;

    vector <double, paracabs::allocator<double, MemTypeDefault>>     vec       (10);
    vector <double, paracabs::allocator<double, MemTypeAccelerator>> vec_accel (10);

    my_copy (vec, vec_accel);

    cout << "Done." << endl;

    return (0);
}