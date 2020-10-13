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
        template <typename type>
        struct Tensor : public Vector<type>
        {
            size_t nrows = 0;
            size_t ncols = 0;
            size_t depth = 0;
            size_t nwarp = ncols;


            ///  Constructor (no argument)
            //////////////////////////////
            inline Tensor ()
            {
                Vector<type>::set_dat ();
            }

            ///  Copy constructor (shallow copy)
            ////////////////////////////////////
            inline Tensor (const Tensor& t)
            {
                nrows = t.nrows;
                ncols = t.ncols;
                depth = t.depth;
                nwarp = t.nwarp;

                Vector<type>::ptr            = t.ptr;
                Vector<type>::allocated      = false;
                Vector<type>::allocated_size = 0;
                Vector<type>::set_dat ();
            }

            ///  Constructor (double argument)
            //////////////////////////////////
            inline Tensor (const size_t nr, const size_t nc, const size_t nd)
            {
                resize (nr, nc, nd);
            }

            ///  Resizing both the std::vector and the allocated memory
            ///    @param[in] size : new size for std::vector
            ///////////////////////////////////////////////////////////
            inline void resize (const size_t nr, const size_t nc, const size_t nd)
            {
                nrows = nr;
                ncols = nc;
                depth = nd;
                nwarp = ncols;

                Vector<type>::vec.resize (nrows*ncols*depth);
                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }

            ///  Access operators
            accel inline type  operator() (const size_t id_r, const size_t id_c, const size_t id_d) const
            {
                return Vector<type>::dat[id_d + depth*(id_c + nwarp*id_r)];
            }

            accel inline type &operator() (const size_t id_r, const size_t id_c, const size_t id_d)
            {
                return Vector<type>::dat[id_d + depth*(id_c + nwarp*id_r)];
            }

            /// Setters for Python
            inline void set_3D_array (py::array_t<type, py::array::c_style | py::array::forcecast> arr)
            {
                py::buffer_info buf = arr.request();

                if (buf.ndim != 3)
                {
                    throw std::runtime_error("Number of dimensions must be 3.");
                }

                type* buf_ptr = (type*) buf.ptr;

                nrows = buf.shape[0];
                ncols = buf.shape[1];
                depth = buf.shape[2];
                nwarp = ncols;

                Vector<type>::vec.resize (nrows*ncols*depth);

                for (size_t i = 0; i < nrows*ncols*depth; i++)
                {
                    Vector<type>::vec[i] = buf_ptr[i];
                }

                Vector<type>::copy_vec_to_ptr ();
                Vector<type>::set_dat ();
            }
        };
    }
}
