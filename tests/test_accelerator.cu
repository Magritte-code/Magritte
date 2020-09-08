#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
using namespace paracabs::accelerator;
using paracabs::datatypes::vector3d;
using paracabs::datatypes::array1d;
using paracabs::datatypes::MemTypeDefault;
using paracabs::datatypes::MemTypeAccelerator;
using paracabs::datatypes::my_copy;


typedef vector3d<double, MemTypeDefault>     Vector3d;
typedef vector3d<double, MemTypeAccelerator> Vector3d_accel;

typedef array1d <Vector3d,       MemTypeDefault>     arrvec3d;
typedef array1d <Vector3d_accel, MemTypeAccelerator> arrvec3d_accel;

typedef array1d <double, MemTypeDefault>     arrdouble;
typedef array1d <double, MemTypeAccelerator> arrdouble_accel;


__global__ void printKernel (arrdouble_accel* arr)
{
    for (size_t i = 0; i < arr->size; i++)
    {
        printf("test");
    }
}


int main ()
{
    cout << "Paracabs test accelerator."   << endl;
    cout << "Number of GPUs = " << nGPUs() << endl;


    list_accelerators();

    const size_t size = 10;

//    arrvec3d*       arr_cpu = new arrvec3d       (size);
//    arrvec3d_accel* arr_gpu = new arrvec3d_accel (size);

    arrdouble*       arr_cpu = new arrdouble       (size);
//    arrdouble_accel* arr_gpu = new arrdouble_accel (size);


    for (int i = 0; i < size; i++)
    {
        arr_cpu[i] = i + 4;
    }

//    paracabs::datatypes::my_copy <double> (arr_cpu, arr_gpu);

//    printKernel <<<1, 1>>> (arr_gpu);
//    paracabs::accelerator::synchronize();

    delete arr_cpu;
//    delete arr_gpu;

    cout << "Done." << endl;

    return (0);
}