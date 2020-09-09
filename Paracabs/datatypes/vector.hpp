#pragma once


#include <vector>
#include "accelerator/accelerator.hpp"



namespace paracabs
{
    namespace datatypes
    {
        template <typename type>
        class Vector
        {
            public:
                type*  dat = nullptr;          ///< pointer to vector data
                type*  ptr = nullptr;          ///< pointer to vector data
                bool   allocated = false;      ///< true if ptr is malloc'ed
                size_t allocated_size = 0;     ///< array size

                vector<type> vec;              ///< stl::vector of data

                ///  Constructor (no argument)
                //////////////////////////////
                inline Vector ()
                {
                    set_dat ();
                }

                ///  Setter for ptr
                ///////////////////
                inline void set_dat ()
                {
                    if (copyContextAccelerator()) dat = ptr;
                    else                          dat = vec.data();
                }

                ///  Copy constructor (shallow copy)
                ////////////////////////////////////
                inline Vector (const Vector& v)
                {
                    ptr            = v.ptr;
                    allocated      = false;
                    allocated_size = 0;
                    set_dat ();
                }

                ///  Constructor (single argument)
                ///    @param[in] v : std::vector to construct Vector from
                //////////////////////////////////////////////////////////
                inline Vector (const vector<type>& v) : vec(v)
                {
                    copy_vec_to_ptr ();
                    set_dat ();
                }

                ///  Constructor (single argument)
                ///    @param[in] s : size of the array to construct
                ////////////////////////////////////////////////////
                inline Vector (const size_t s) : vec(s)
                {
                    copy_vec_to_ptr ();
                    set_dat ();
                }

                ///  Constructor (two arguments)
                ///    @param[in] s : size of the array to construct
                ///    @param[in] i : initializer values
                ////////////////////////////////////////////////////
                inline Vector (const size_t s, const type i) : vec(s, i)
                {
                    copy_vec_to_ptr ();
                    set_dat ();
                }

                ///  Destructor
                ///////////////
                inline ~Vector () {deallocate();}

                ///  Memory allocator
                ///    @param[in] size : number of elements
                ///////////////////////////////////////////
                inline void allocate (const size_t size)
                {
                    #if PARACABS_USE_ACCELERATOR
                        if (allocated_size != size)
                        {
                            if (allocated) deallocate();
                            ptr = (type*) paracabs::accelerator::malloc (size*sizeof(type));
                            allocated = true;
                            allocated_size = size;
                            set_dat ();
                        }
                    #endif
                }

                ///  Memory deallocator
                ///////////////////////
                inline void deallocate ()
                {
                    #if PARACABS_USE_ACCELERATOR
                        if (allocated)
                        {
                            paracabs::accelerator::free (ptr);
                            allocated = false;
                            allocated_size = 0;
                        }
                    #endif
                }

                ///  Copier for data from std::vector to allocated memory
                /////////////////////////////////////////////////////////
                inline void copy_vec_to_ptr ()
                {
                    #if PARACABS_USE_ACCELERATOR
                        allocate (vec.size());
                        paracabs::accelerator::memcpy_to_accelerator (ptr, vec.data(), vec.size()*sizeof(type));
                        set_dat ();
                    #endif
                }

                ///  Copier for data from allocated memory tp std::vector
                /////////////////////////////////////////////////////////
                inline void copy_ptr_to_vec ()
                {
                    #if PARACABS_USE_ACCELERATOR
                        vec.resize (allocated_size);
                        paracabs::accelerator::memcpy_from_accelerator (vec.data(), ptr, vec.size()*sizeof(type));
                        set_dat ();
                    #endif
                }

                ///  Resizing both the std::vector and the allocated memory
                ///    @param[in] size : new size for std::vector
                ///////////////////////////////////////////////////////////
                inline void resize (const size_t size)
                {
                    vec.resize (size);
                    copy_vec_to_ptr ();
                    set_dat ();
                }

                ///  Access operators
                accel inline type  operator[] (const size_t id) const {return dat[id];}
                accel inline type &operator[] (const size_t id)       {return dat[id];}
        };


//        template <class C>
//        struct MM : public C
//        {
//            MM* ptr = nullptr;
//
//            MM ()
//            {
//                ptr = (MM*) paracabs::accelerator::malloc (sizeof(MM));
//
//                paracabs::accelerator::memcpy_to_accelerator (ptr, this, sizeof(MM));
//            }
//
//            MM (const MM<C>& obj)
//            {
//                this = ptr
//            }
//        };

    }
}
