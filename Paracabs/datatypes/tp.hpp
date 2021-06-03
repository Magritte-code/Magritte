#pragma once


#include "datatypes/vector.hpp"
#include "accelerator/accelerator.hpp"
#include "multi_threading/multi_threading.hpp"


namespace paracabs
{
    namespace datatypes
    {
        /// TP: a thread private 0-index data structure
        ///////////////////////////////////////////////
        template <typename type, typename XThreads>
        struct TP : public Vector<type>, XThreads
        {
            ///  Constructor (no argument)
            //////////////////////////////
            inline TP ()
            {
                Vector<type>::vec.resize (XThreads::tot_nthreads());
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Copy constructor (shallow copy)
            ////////////////////////////////////
            inline TP (const TP& v)
            {
                Vector<type>::ptr            = v.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Access operators
            accel inline type  operator() () const
            {
                return Vector<type>::dat[XThreads::thread_id()];
            }

            accel inline type &operator() ()
            {
                return Vector<type>::dat[XThreads::thread_id()];
            }

            accel inline type  operator() (const size_t t) const
            {
                return Vector<type>::dat[t];
            }

            accel inline type &operator() (const size_t t)
            {
                return Vector<type>::dat[t];
            }
        };
    }
}
