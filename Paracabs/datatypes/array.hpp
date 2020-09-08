#pragma once


#include "accelerator/accelerator.hpp"


namespace paracabs
{
    namespace datatypes
    {
        template <typename type, class MemType>
        class Array : public MemType
        {
            private:
                type* data;        ///< array data
//                bool  allocated;   ///< true if memory was already allocated?

            public:
                const size_t size;   ///< array size

                ///  Constructor (single argument)
                ///    @param[in] s : size of the array to construct
                ////////////////////////////////////////////////////
                accel inline Array (const size_t s) : size (s)
                {
                    data = (type*) this->malloc (size*sizeof(type));
                }


//                Array ()
//                {
//                     // This only works if Array has MemTypeDefault

//                    data = (type*) this->malloc (size*sizeof(type));
//
//                    for (size_t n = 0; n < size; n++)
//                    {
//                        data[n] = vec[n];
//                    }
//                }

//                stl::vector<type> to_stl_vector () const
//                {
//                    vector<type> vec (size);
//                }

                ///  Destructor
                ///////////////
                accel inline ~Array ()
                {
                    this->free (data);
                }

                ///  Access operators
                accel inline type  operator[] (const size_t i) const {return data[i];}
                accel inline type &operator[] (const size_t i)       {return data[i];}

//                accel inline bool is_malloced () {return allocated;}
        };
    }
}