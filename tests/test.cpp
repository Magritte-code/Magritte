#include <iostream>
using std::cout;
using std::endl;

#include "paracabs.hpp"
using namespace paracabs;

struct Model_Vector
{
    const size_t n = 10;

    datatypes::Vector<double> a;
    datatypes::Vector<double> b;
    datatypes::Vector<double> c;

    int fun ()
    {
        accelerator::nblocks  = 1;
        accelerator::nthreads = 1;
        
        a.resize(n);
        b.resize(n);
        c.resize(n);

        for (size_t i = 0; i < n; i++)
        {
            a[i] = 0.5;
            b[i] = 1.5;
            c[i] = 0.0;

            // printf("a[i] = %lf\n", a[i]);
            // printf("b[i] = %lf\n", b[i]);
            // printf("c[i] = %lf\n", c[i]);
        }

        a.copy_vec_to_ptr();
        b.copy_vec_to_ptr();

        accelerated_for (i, n,
        {
            c[i] = a[i] + b[i];

            // printf("a[i] = %lf\n", a[i]);
            // printf("b[i] = %lf\n", b[i]);
            // printf("c[i] = %lf\n", c[i]);
        })

        accelerator::synchronize();

        c.copy_ptr_to_vec();

        for (size_t i = 0; i < n; i++)
        {
            cout << i << "  " << c[i] << endl;

            if (c[i] != 2.0)
            {
                // throw std::runtime_error("ERROR: c[i] != 2.0");
            }
        }

        return (0);
    }
};


struct Model_VectorTP
{
    const size_t n = 10;

    datatypes::VectorTP<double, accelerator::AcceleratorThreads> a;
    datatypes::VectorTP<double, accelerator::AcceleratorThreads> b;
    datatypes::VectorTP<double, accelerator::AcceleratorThreads> c;

    int fun ()
    {
        accelerator::nthreads = 1;
        accelerator::nblocks  = 2;

        a.resize(n);
        b.resize(n);
        c.resize(n);

        cout << "a.tot_nthreads() = " << a.tot_nthreads() << endl;
        cout << "a.size()         = " << a.size()         << endl;

        for (size_t t = 0; t < a.tot_nthreads(); t++)
        {
            for (size_t i = 0; i < n; i++)
            {
                a(t,i) = 0.5;
                b(t,i) = 1.5;
                c(t,i) = 7.0;
            }
        }

        a.copy_vec_to_ptr();
        b.copy_vec_to_ptr();


        accelerated_for (i, n,
        {
            // printf("id = %lf", paracabs::multi_threading::thread_id());
            // printf("test\n");
            printf("id = %ld\n", a.thread_id());

            c[i] = a[i] + b[i];

            // printf("a[i] = %lf\n", a[i]);
            // printf("b[i] = %lf\n", b[i]);
            // printf("c[i] = %lf\n", c[i]);
        })


        accelerator::synchronize();

        c.copy_ptr_to_vec();


        for (size_t t = 0; t < a.tot_nthreads(); t++)
        {
            for (size_t i = 0; i < n; i++)
            {
                cout << i << "  " << c(t,i) << endl;

                if (c(t,i) != 2.0)
                {
                    //throw std::runtime_error("ERROR: c[i] != 2.0");
                }
            }
        }

        return (0);
    }
};


int test_Vector_accelerated_for ()
{
    accelerator::list_accelerators();

    Model_Vector model;
    model.fun();

    return (0);
}


int test_VectorTP_accelerated_for ()
{
    accelerator::list_accelerators();

    Model_VectorTP model;
    model.fun();

    return (0);
}


int test_accelerated_for_outside_class ()
{
    const size_t n = 100;

    accelerator::nblocks  = 1;
    accelerator::nthreads = 1;

    datatypes::Vector<double> a (n, 0.5);
    datatypes::Vector<double> b (n, 1.5);
    datatypes::Vector<double> c (n, 0.0);

    accelerator::list_accelerators();

    accelerated_for_outside_class (i, n,
    {
        c[i] = a[i] + b[i];
    })

    accelerator::synchronize();

    c.copy_ptr_to_vec();

    for (size_t i = 0; i < n; i++)
    {
        if (c[i] != 2.0)
        {
            throw std::runtime_error("ERROR: c[i] != 2.0");
        }
    }

    return (0);
}





int main ()
{
    //test_Vector_accelerated_for();

    test_VectorTP_accelerated_for();

    // test_accelerated_for_outside_class();

    // multi_threading::set_n_threads_avail(2);
    // cout << multi_threading::n_threads_avail() << endl;

    return (0);
}
