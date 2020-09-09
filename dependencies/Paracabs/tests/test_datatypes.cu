#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
namespace pc = paracabs::datatypes;

using Vector3D = pc::Vector3D<double>;

template <typename type>
using Array     = pc::Array <type, pc::MemTypeDefault>;
template <typename type>
using Array_acc = pc::Array <type, pc::MemTypeAccelerator>;



struct Model
{
    size_t length;
    pc::Vector <double> d;

    Model (const size_t s) : length(s), d(s, 8.44) {}

//    Model (const Model& m)
//    {
//        length = m.length;
//        d      = pc::Vector<double> (m.d);
//    }

    inline void f ()
    {
        for (int i = 0; i < length; i++)
        {
            d[i] = 2.29;
        }
    }

    inline void p ()
    {
        for (int i = 0; i < length; i++)
        {
            cout << d[i] << endl;
        }
    }
};


struct Test
{
    size_t length;
    pc::Vector <double> a;

    Model model;


    Test (const size_t s) : length (s), model(3) {a.resize(s);}

//    double* ptr;
//    double* ptr_host;
//    double* ptr_device;
//
//    accel inline type  operator[] (const size_t id) const {return ptr[id];}
//    accel inline type &operator[] (const size_t id)       {return ptr[id];}
//
//    Test (const Test& test)
//    {
//        if (accelerator_code) ptr = ptr_host;
//        else                  ptr = ptr_device;
//    }


    void test ()
    {
        for (int i = 0; i < length; i++)
        {
            a[i] = 3.14;
            cout << a[i] << endl;
        }

        a.copy_vec_to_ptr ();

        accelerated_for (i, length, 1, 1,
        {
            a[i] -= 3;
        })

        a.copy_ptr_to_vec ();

        cout << endl;

        for (int i = 0; i < length; i++)
        {
            cout << a[i] << endl;
        }
    }

    void test2 ()
    {
        model.p();

//        a.copy_vec_to_ptr ();

        accelerated_for (i, model.length, 1, 1,
        {
            model.d[i] -= 3;
        })

        model.d.copy_ptr_to_vec ();

        model.p();
    }

};


//__global__ void addKernel ()
//{
//
//}


int main ()
{
    cout << "Paracabs test datatypes." << endl;

//    Vector3D v1 (1.0, 2.0, 3.0);
//    Vector3D v2 (4.0, 5.0, 6.0);
//    Vector3D v3 = v1 + v2;
//
//    v1.print();
//    v2.print();
//    v3.print();
//
//    cout << v1.dot(v2) << endl;
//
//
//    v1 += v2;
//
//    v1.print();
//    v2.print();
//
//
//    Vector3D v4 = 3.14;
//
//    v4.print();
//
//    (v4 + 1).print();
//
//    v4.print();
//
//    v4 = 7.12;
//
//    v4.print();
//
//
//    const size_t size = 10;

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

    Test t (10);

//    t.test();

//    size_t length = 15;
//    pc::Vector <double> a (length);
//
//    for (int i = 0; i < length; i++)
//    {
//        a[i] = 3.14;
//        cout << a[i] << endl;
//    }
//
//    a.copy_vec_to_ptr ();
//
//    accelerated_for_outside_class (i, length, 1, 1,
//    {
//        a[i] -= 3;
//    })
//
//    a.copy_ptr_to_vec ();
//
//    cout << endl;
//
//    for (int i = 0; i < length; i++)
//    {
//        cout << a[i] << endl;
//    }

    cout << "Done." << endl;

//    t.model.p();
//    t.model.f();
//    t.model.p();

    t.test2();



    return (0);
}