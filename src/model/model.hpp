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


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;
    Radiation      radiation;

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
    int compute_spectral_discretisation           (const Real width);
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

    Double1 error_max;
    Double1 error_mean;

    pc::multi_threading::ThreadPrivate<Vector<Real>> a;
    pc::multi_threading::ThreadPrivate<Vector<Real>> b;
    pc::multi_threading::ThreadPrivate<Vector<Real>> c;


    int set()
    {
       for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
       {
           const Size N = 1000000;

           cout << "Setting thread " << i << endl;

           a(i).resize(N);
           b(i).resize(N);
           c(i).resize(N);

           for (Size n = 0; n < N; n++)
           {
               a(i)[n] = 1.0;
               b(i)[n] = 2.0;
               c(i)[n] = 5.0;
           }
       }

       return (0);
    }

    Vector<Real> add ()
    {
        cout << "c.size() = " << c().vec.size() << endl;

        accelerated_for (i, c().vec.size(), nblocks, nthreads,
        {
            c()[i] = a()[i] + b()[i];

            // cout << "c[" << i << "] = " << c()[i] << endl;
        })

        return c().vec;
    }


};
