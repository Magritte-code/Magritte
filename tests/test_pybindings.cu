#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
namespace pc = paracabs::datatypes;

using Vector3D = pc::Vector3D<double, pc::MemTypeDefault>;

template <typename type>
using Array     = pc::Array <type, pc::MemTypeDefault>;
template <typename type>
using Array_acc = pc::Array <type, pc::MemTypeAccelerator>;


struct test
{
    Array     <double>*   nums;
    Array_acc <double>*   nums_acc;
    Array     <Vector3D>* vecs;
    Array_acc <Vector3D>* vecs_acc;


};


__global__ void addKernel ()
{

}


int main ()
{
    cout << "Paracabs test datatypes." << endl;

    Vector3D v1 (1.0, 2.0, 3.0);
    Vector3D v2 (4.0, 5.0, 6.0);
    Vector3D v3 = v1 + v2;

    v1.print();
    v2.print();
    v3.print();

    cout << v1.dot(v2) << endl;


    v1 += v2;

    v1.print();
    v2.print();


    Vector3D v4 = 3.14;

    v4.print();

    (v4 + 1).print();

    v4.print();

    v4 = 7.12;

    v4.print();


    const size_t size = 10;

//    array1d <vector3d <double>, MemTypeDefault>     arr1 (size);
//    array1d <vector3d <double>, MemTypeAccelerator> arr2 (size);
//
//    for (size_t i = 0; i < size; i++)
//    {
//        arr1[i] = 1.0;
//    }
//
//    for (size_t i = 0; i < size; i++)
//    {
//        arr1[i].print();
//    }
//
//    for (size_t i = 0; i < size; i++)
//    {
//        arr2[i] = 1.0;
//    }
//
//    for (size_t i = 0; i < size; i++)
//    {
//        arr2[i].print();
//    }


    cout << "Done." << endl;

    return (0);
}