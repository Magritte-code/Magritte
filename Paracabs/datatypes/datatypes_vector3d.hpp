#pragma once


namespace datatypes
{
    template <typename type>
    struct vector3d
    {
        type data[3];


        vector3d (const type x, const type y, const type z)
        {
            data[0] = x;
            data[1] = y;
            data[2] = z;
        }


        inline vector3d<type> operator+ (const vector3d<type> vec) const
        {
            const type x = data[0] + vec.data[0];
            const type y = data[1] + vec.data[1];
            const type z = data[2] + vec.data[2];

            return vector3d<type> (x, y, z);
        }


        inline vector3d<type> operator- (const vector3d<type> vec) const
        {
            const type x = data[0] - vec.data[0];
            const type y = data[1] - vec.data[1];
            const type z = data[2] - vec.data[2];

            return vector3d<type> (x, y, z);
        }


        inline void operator+= (const vector3d<type> vec)
        {
            data[0] += vec.data[0];
            data[1] += vec.data[1];
            data[2] += vec.data[2];
        }


        inline void operator-= (const vector3d<type> vec)
        {
            data[0] -= vec.data[0];
            data[1] -= vec.data[1];
            data[2] -= vec.data[2];
        }


        inline type x() const {return data[0];}
        inline type y() const {return data[1];}
        inline type z() const {return data[2];}


        inline void print () const {printf("%lf, %lf, %lf\n", x(), y(), z());}
    };


    template <typename type>
    inline type dot (const vector3d<type> &a, const vector3d<type> &b)
    {
        return   a.data[0] * b.data[0]
               + a.data[1] * b.data[1]
               + a.data[2] * b.data[2];
    }

}

