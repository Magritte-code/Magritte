#include <iostream>
using std::cout;
using std::endl;

#include "paracabs.hpp"
using namespace paracabs;

struct Model
{
    const size_t n = 10;

    datatypes::Vector<double> a;
    datatypes::Vector<double> b;
    datatypes::Vector<double> c;

    int fun ()
    {
        a.resize(n);
        b.resize(n);
        c.resize(n);

        for (size_t i = 0; i < n; i++)
        {
            a[i] = 0.5;
            b[i] = 1.5;
            c[i] = 0.0;

            printf("a[i] = %lf\n", a[i]);
            printf("b[i] = %lf\n", b[i]);
            printf("c[i] = %lf\n", c[i]);
        }

        a.copy_vec_to_ptr();
        b.copy_vec_to_ptr();

        const size_t nblocks  = 1;
        const size_t nthreads = 1;
        
        accelerated_for (i, n, nblocks, nthreads,
        {
            c[i] = a[i] + b[i];

            printf("a[i] = %lf\n", a[i]);
            printf("b[i] = %lf\n", b[i]);
            printf("c[i] = %lf\n", c[i]);
        })

        accelerator::synchronize();

        c.copy_ptr_to_vec();

        for (size_t i = 0; i < n; i++)
        {
            cout << i << "  " << c[i] << endl;

            if (c[i] != 2.0)
            {
                //throw std::runtime_error("ERROR: c[i] != 2.0");
            }
        }

        return (0);
    }
};


int test_accelerated_for ()
{
    accelerator::list_accelerators();

    Model model;
    model.fun();

    cout << "Done" << endl;
    
    return (0);
}


int test_accelerated_for_outside_class ()
{
    const size_t n = 100;

    const size_t nblocks  = 1;
    const size_t nthreads = 1;

    datatypes::Vector<double> a (n, 0.5);
    datatypes::Vector<double> b (n, 1.5);
    datatypes::Vector<double> c (n, 0.0);

    accelerator::list_accelerators();

    accelerated_for_outside_class (i, n, nblocks, nthreads,
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

    cout << "Done" << endl;
    
    return (0);
}





int main ()
{
    test_accelerated_for();

    test_accelerated_for_outside_class();

    return (0);
}
