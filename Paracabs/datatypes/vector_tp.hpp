#pragma once


#include "datatypes/vector.hpp"
#include "accelerator/accelerator.hpp"
#include "multi_threading/multi_threading.hpp"


namespace paracabs
{
    namespace datatypes
    {
        /// VectorTP: a thread private 1-index data structure
        /////////////////////////////////////////////////////
        template <typename type>
        struct VectorTP : public Vector<type>
        {
            size_t size     = 0;
            size_t nthreads = paracabs::multi_threading::n_threads_avail();


            ///  Constructor (no argument)
            //////////////////////////////
            inline VectorTP ()
            {
                Vector<type>::set_dat ();
            }

            ///  Copy constructor (shallow copy)
            ////////////////////////////////////
            inline VectorTP (const VectorTP& v)
            {
                size     = v.size;
                nthreads = v.nthreads;

                Vector<type>::ptr            = v.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Constructor (double argument)
            //////////////////////////////////
            inline VectorTP (const size_t s)
            {
                VectorTP<type>::resize (s);
            }

            ///  Resizing both the std::vector and the allocated memory
            ///    @param[in] size : new size for std::vector
            ///////////////////////////////////////////////////////////
            inline void resize (const size_t s)
            {
                size = s;

                Vector<type>::vec.resize (size*nthreads);
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Access operators
            accel inline type  operator[] (const size_t id) const
            {
                return Vector<type>::dat[id + size*paracabs::multi_threading::thread_id()];
            }

            accel inline type &operator[] (const size_t id)
            {
                return Vector<type>::dat[id + size*paracabs::multi_threading::thread_id()];
            }

            accel inline type  operator() (const size_t t, const size_t id) const
            {
                return Vector<type>::dat[id + size*t];
            }

            accel inline type &operator() (const size_t t, const size_t id)
            {
                return Vector<type>::dat[id + size*t];
            }
        };
    }
}
