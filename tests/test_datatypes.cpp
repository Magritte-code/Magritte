#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
using paracabs::datatypes::vector3d;
using paracabs::datatypes::array1d;


int main ()
{
    cout << "Paracabs test datatypes." << endl;

    vector3d <double> v1 (1.0, 2.0 , 3.0);
    vector3d <double> v2 (4.0, 5.0 , 6.0);

    vector3d <double> v3 = v1 + v2;

    v1.print();
    v2.print();
    v3.print();

    cout << dot(v1,v2) << endl;


    v1 += v2;

    v1.print();
    v2.print();

    const size_t size = 10;

    array1d <vector3d <double>, paracabs::datatypes  ::MemoryManagement> arr1 (size);
    array1d <vector3d <double>, paracabs::accelerator::MemoryManagement> arr2 (size);




    cout << "Done." << endl;

    return (0);
}