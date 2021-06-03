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
        template <typename type, typename XThreads>
        struct VectorTP : public Vector<type>, XThreads
        {
            size_t vec_size = 0;


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
                vec_size = v.vec_size;

                Vector<type>::ptr            = v.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Constructor (double argument)
            //////////////////////////////////
            inline VectorTP (const size_t s)
            {
                VectorTP<type, XThreads>::resize (s);
            }

            ///  Resizing both the std::vector and the allocated memory
            ///    @param[in] size : new size for std::vector
            ///////////////////////////////////////////////////////////
            inline void resize (const size_t s)
            {
                vec_size = s;

                Vector<type>::vec.resize (vec_size*XThreads::tot_nthreads());
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Indexing
            /////////////
            accel inline size_t index (const size_t t, const size_t id) const
            {
                return id + vec_size*t;
            }

            ///  Access operators
            accel inline type  operator[] (const size_t id) const
            {
                return Vector<type>::dat[index(XThreads::thread_id(), id)];
            }

            accel inline type &operator[] (const size_t id)
            {
                return Vector<type>::dat[index(XThreads::thread_id(), id)];
            }

            accel inline type  operator() (const size_t t, const size_t id) const
            {
                return Vector<type>::dat[index(t, id)];
            }

            accel inline type &operator() (const size_t t, const size_t id)
            {
                return Vector<type>::dat[index(t, id)];
            }
        };
    }
}
