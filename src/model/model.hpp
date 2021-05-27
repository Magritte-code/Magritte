#pragma once


#include "io/io.hpp"
#include "io/python/io_python.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"
#include "radiation/radiation.hpp"
#include "image/image.hpp"
#include "paracabs.hpp"


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;
    Radiation      radiation;
    vector<Image>  images;

    enum SpectralDiscretisation {None, SD_Lines, SD_Image}
         spectralDiscretisation = None;

    Model () {};
    Model (const string name)
    {
        parameters.model_name() = name;
        read ();
    }

    void read  (const Io& io);
    void write (const Io& io) const;

    void read  ()       {read  (IoPython ("hdf5", parameters.model_name()));};
    void write () const {write (IoPython ("hdf5", parameters.model_name()));};

    int compute_inverse_line_widths               ();
    int compute_spectral_discretisation           ();
    int compute_spectral_discretisation           (
        const Real width );
    int compute_spectral_discretisation           (
        const long double nu_min,
        const long double nu_max );
    int compute_LTE_level_populations             ();
    int compute_radiation_field                   ();
    int compute_radiation_field_feautrier_order_2 ();
    int compute_radiation_field_shortchar_order_0 ();
    int compute_Jeff                              ();
    int compute_level_populations_from_stateq     ();
    int compute_level_populations                 (
        // const Io   &io,
        const bool  use_Ng_acceleration,
        const long  max_niterations     );
    int compute_image                             (const Size ray_nr);

    Double1 error_max;
    Double1 error_mean;

    pc::multi_threading::ThreadPrivate<Vector<Real>> a;
    pc::multi_threading::ThreadPrivate<Vector<Real>> b;
    pc::multi_threading::ThreadPrivate<Vector<Real>> c;

    Vector<Real> x;
    Vector<Real> y;
    Vector<Real> z;


    int set()
    {
        //for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
        //{
        //    const Size N = 100;

        //    cout << "Setting thread " << i << endl;

        //    a(i).resize(N);
        //    b(i).resize(N);
        //    c(i).resize(N);

        //    for (Size n = 0; n < N; n++)
        //    {
        //        a(i)[n] = 1.0;
        //        b(i)[n] = 2.0;
        //        c(i)[n] = 5.0;
        //    }
        //}

        const Size N = 100;

        x.resize(N);
        y.resize(N);
        z.resize(N);

        for (Size n = 0; n < N; n++)
        {
            x[n] = 1.0;
            y[n] = 2.0;
            z[n] = 5.0;
        }

        x.copy_vec_to_ptr(); 
        y.copy_vec_to_ptr(); 
        z.copy_vec_to_ptr(); 

        return (0);
    }



    Vector<Real> add ()
    {
        cout << "c.size() = " << c().vec.size() << endl;

        // accelerated_for (i, 100, 1, 1,
        // {
        //     c()[i] = a()[i] + b()[i];
        // })

        copyContextAccelerator() = true;                                          
        auto lambda = [=, *this] __device__ (size_t i) mutable                    
        {		                                                                  
            // c()[i] = a()[i] + b()[i];
            x[i] = y[i] + z[i];
        };							

        decltype(lambda)* lambda_ptr =                                            
            (decltype(lambda)*) paracabs::accelerator::malloc (sizeof(lambda));   
        
        paracabs::accelerator::memcpy_to_accelerator                              
            (lambda_ptr, &lambda, sizeof(lambda));                                
        
        apply_lambda <<<1, 1>>> (100, lambda_ptr);                 
        copyContextAccelerator() = false;                                         

        pc::accelerator::synchronize();

        c().copy_ptr_to_vec();

        for (int i = 0; i < 100; i++)
        {
            // cout << "c["<< i <<"] = " << c()[i] << endl;
            cout << "x["<< i <<"] = " << x[i] << endl;
        }

        //return c().vec;
        return x.vec;
    }

    // Kernel approach
    Matrix<Real> eta;
    Matrix<Real> chi;

    Matrix<Real> boundary_condition;

    int set_eta_and_chi       ();
    int set_boundary_condition();
};
