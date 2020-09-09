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

    inline void print ()
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


    Test (const size_t s): length(s), model(3) {a.resize(s);}

    void test ()
    {
        for (size_t i = 0; i < length; i++)
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

        for (size_t i = 0; i < length; i++)
        {
            cout << a[i] << endl;
        }
    }

    void test2 ()
    {
        model.print();

        accelerated_for (i, model.length, 1, 1,
        {
//            cout << model.d[i] << "  ";
            model.d[i] -= 3;
//            cout << model.d[i] << endl;
        })

        model.d.copy_ptr_to_vec ();

        model.print();
    }

};



int main ()
{
    cout << "Paracabs test datatypes." << endl;

    Test t (10);
    t.test2();

    cout << "Done." << endl;


    return (0);
}