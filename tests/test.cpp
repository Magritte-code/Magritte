#include <iostream>
using std::cout;
using std::endl;

#include "paracabs.hpp"


int main ()
{
    cout << "Done." << endl;

    double arr[66];

    cout << sizeof(arr) / sizeof(double) << endl;

    class Vec
    {
        double arr1[5];
        double arr2[2];
        double scl;
    } vec;

    cout << sizeof(vec) / sizeof(double) << endl;

    return (0);
}