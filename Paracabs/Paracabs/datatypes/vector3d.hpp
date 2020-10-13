#pragma once

#include "accelerator/accelerator.hpp"


namespace paracabs
{
    namespace datatypes
    {
        template <typename type>
        class Vector3D
        {
            public:
                type data[3];   ///< vector components

                ///  Constructor
                /////////////////////////////////
                accel inline Vector3D () {}

                ///  Constructor
                ///    @param[in] x : x component
                ///    @param[in] y : y component
                ///    @param[in] z : z component
                /////////////////////////////////
                accel inline Vector3D (const type x, const type y, const type z)
                {
                    data[0] = x;
                    data[1] = y;
                    data[2] = z;
                }

                ///  Constructor (single argument)
                ///    @param[in] v : x, y and z component
                //////////////////////////////////////////
                accel inline Vector3D (const type v)
                {
                    *this = Vector3D (v, v, v);
                }

//        inline vector3d &operator= (const vector3d &&rhs)
//        {
//            data = rhs.data;
//            return *this;
//        };

//        inline vector3d &operator= (const vector3d &rhs)
//        {
//            data = rhs.data;
//            return *this;
//        };

//        vector3d () = default;
//        vector3d (const vector3d &rhs ) : data (rhs.data){};  // compiles in movaps
//        vector3d (const vector3d &&rhs) : data (rhs.data){};


//        inline void operator= (const vector3d& vec)
//        {
//            &data = &vec.data;
//        }


                ///  assignment operator
                ////////////////////////
                accel inline void operator= (const type v)
                {
                    *this = Vector3D (v);
                }

                ///  Addition operator
                //////////////////////
                accel inline Vector3D operator+ (const Vector3D& vec) const
                {
                    const type x = data[0] + vec.data[0];
                    const type y = data[1] + vec.data[1];
                    const type z = data[2] + vec.data[2];

                    return Vector3D (x, y, z);
                }

                ///  Subtraction operator
                /////////////////////////
                accel inline Vector3D operator- (const Vector3D& vec) const
                {
                    const type x = data[0] - vec.data[0];
                    const type y = data[1] - vec.data[1];
                    const type z = data[2] - vec.data[2];

                    return Vector3D (x, y, z);
                }

                ///  Add assignment operator
                ////////////////////////////
                accel inline void operator+= (const Vector3D& vec)
                {
                    data[0] += vec.data[0];
                    data[1] += vec.data[1];
                    data[2] += vec.data[2];
                }

                ///  Subtraction assignment operator
                ////////////////////////////////////
                accel inline void operator-= (const Vector3D& vec)
                {
                    data[0] -= vec.data[0];
                    data[1] -= vec.data[1];
                    data[2] -= vec.data[2];
                }

                ///  (Euclidean) dot product
                ///    @param[in] vec : vector to take dot product with
                ///    @returns the dot product between this vector and vec
                ///////////////////////////////////////////////////////////
                accel inline type dot (const Vector3D& vec) const
                {
                    return   data[0] * vec.data[0]
                           + data[1] * vec.data[1]
                           + data[2] * vec.data[2];
                }

                ///  Squared norm (dot product with itself)
                ///////////////////////////////////////////
                accel inline type squaredNorm () const
                {
                    return dot (*this);
                }


                accel inline type x () const {return data[0];}
                accel inline type y () const {return data[1];}
                accel inline type z () const {return data[2];}

                accel inline void print () const {printf ("%le, %le, %le\n", x(), y(), z());}
        };

    }
}
