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


namespace sycl = cl::sycl;


struct Model
{
    size_t length;
    pc::Vector <double> d;
    //double* d;

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
    size_t length = 3;
    
    //pc::Vector <double> a = pc::Vector<double>(3);
    //double* a;

    Model model;


    Test (const size_t s): length(s), model(3) {};

    //Test (const Test& test):
    //    length(test.length),
    //    model(test.model),
    //    a (test.a) {};



    //void test ()
    //{
    //    for (size_t i = 0; i < length; i++)
    //    {
    //        a[i] = 3.14;
    //        cout << a[i] << endl;
    //    }

    //    a.copy_vec_to_ptr ();

    //    accelerated_for (i, length, 1, 1,
    //    {
    //        a[i] -= 3;
    //    })

    //    a.copy_ptr_to_vec ();

    //    cout << endl;

    //    for (size_t i = 0; i < length; i++)
    //    {
    //        cout << a[i] << endl;
    //    }
    //}

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


sycl::default_selector device_selector;
    
sycl::queue paracabs::accelerator::acceleratorQueue = sycl::queue (device_selector);


int main ()
{
    cout << "Paracabs test datatypes." << endl;

    paracabs::accelerator::list_accelerators();


    

    Test t(3);
    t.test2();

    //sycl::float4 a = { 1.0, 2.0, 3.0, 4.0 };
    //sycl::float4 b = { 4.0, 3.0, 2.0, 1.0 };
    //sycl::float4 c = { 0.0, 0.0, 0.0, 0.0 };

    //sycl::default_selector device_selector;Y

    //sycl::queue queue(device_selector);

    //std::cout << "Running on "
    //          << queue.get_device().get_info<sycl::info::device::name>()
    //          << "\n";

    //{ // start of scope, ensures data coped back to host
    //    sycl::buffer<sycl::float4, 1> a_sycl(&a, sycl::range<1>(1));
    //    sycl::buffer<sycl::float4, 1> b_sycl(&b, sycl::range<1>(1));
    //    sycl::buffer<sycl::float4, 1> c_sycl(&c, sycl::range<1>(1));

    //    queue.submit([&] (sycl::handler& cgh) {
    //        auto a_acc = a_sycl.get_access<sycl::access::mode::read>(cgh);
    //        auto b_acc = b_sycl.get_access<sycl::access::mode::read>(cgh);
    //        auto c_acc = c_sycl.get_access<sycl::access::mode::discard_write>(cgh);

    //        cgh.single_task<class vector_addition>([=] () {
    //            c_acc[0] = a_acc[0] + b_acc[0];
    //        });
    //    });

    //} // end of scope, ensures data copied back to host


    cout << "Done." << endl;


    return (0);
}
