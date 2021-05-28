#pragma once


#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
namespace py = pybind11;


#include "datatypes/vector.hpp"
#include "accelerator/accelerator.hpp"


namespace paracabs
{
    namespace datatypes
    {
        /// MatrixTP: a thread private 2-index data structure
        /////////////////////////////////////////////////////
        template <typename type>
        struct MatrixTP : public Vector<type>
        {
            size_t nrows = 0;
            size_t ncols = 0;
            size_t nthreads = paracabs::multi_threading::n_threads_avail();


            ///  Constructor (no argument)
            //////////////////////////////
            inline MatrixTP ()
            {
                Vector<type>::set_dat ();
            }

            ///  Copy constructor (shallow copy)
            ////////////////////////////////////
            inline MatrixTP (const MatrixTP& m)
            {
                nrows    = m.nrows;
                ncols    = m.ncols;
                nthreads = m.nthreads;

                Vector<type>::ptr            = m.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Constructor (double argument)
            //////////////////////////////////
            inline Tensor (const size_t nr, const size_t nc)
            {
                MatrixTP::resize (nr, nc);
            }

            ///  Resizing both the std::vector and the allocated memory
            ///    @param[in] size : new size for std::vector
            ///////////////////////////////////////////////////////////
            inline void resize (const size_t nr, const size_t nc)
            {
                nrows = nr;
                ncols = nc;

                Vector<type>::vec.resize (nrows*ncols*nthreads);
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Access operators
            accel inline type  operator() (const size_t id_r, const size_t id_c) const
            {
                return Vector<type>::dat[id_c + ncols*(id_r + nrows*paracabs::multi_threading::thread_id())];
            }

            accel inline type &operator() (const size_t id_r, const size_t id_c)
            {
                return Vector<type>::dat[id_c + ncols*(id_r + nrows*paracabs::multi_threading::thread_id())];
            }

            accel inline type  operator() (const size_t t, const size_t id_r, const size_t id_c) const
            {
                return Vector<type>::dat[id_c + ncols*(id_r + nrows*t)];
            }

            accel inline type &operator() (const size_t t, const size_t id_r, const size_t id_c)
            {
                return Vector<type>::dat[id_c + ncols*(id_r + nrows*t)];
            }
        };
    }
}
