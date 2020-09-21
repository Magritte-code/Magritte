#pragma once


#include "datatypes/vector.hpp"
#include "accelerator/accelerator.hpp"



namespace paracabs
{
    namespace datatypes
    {
        template <typename type>
        struct Matrix : public Vector<type>
        {
            size_t nrows = 0;
            size_t ncols = 0;

            size_t nwarp = ncols;


            ///  Constructor (no argument)
            //////////////////////////////
            inline Matrix ()
            {
                Vector<type>::set_dat ();
            }

            ///  Copy constructor (shallow copy)
            ////////////////////////////////////
            inline Matrix (const Matrix& m)
            {
                nrows = m.nrows;
                ncols = m.ncols;
                nwarp = m.nwarp;

                Vector<type>::ptr            = m.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Resizing both the std::vector and the allocated memory
            ///    @param[in] size : new size for std::vector
            ///////////////////////////////////////////////////////////
            inline void resize (const size_t nr, const size_t nc)
            {
                nrows = nr;
                ncols = nc;

                nwarp = ncols;

                Vector<type>::vec.resize (nrows*ncols);
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Access operators
            accel inline type  operator() (const size_t id_r, const size_t id_c) const
            {
                return Vector<type>::dat[id_c + nwarp*id_r];
            }

            accel inline type &operator() (const size_t id_r, const size_t id_c)
            {
                return Vector<type>::dat[id_c + nwarp*id_r];
            }
        };
    }
}
