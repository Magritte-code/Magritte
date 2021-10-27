#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Rays
{
    Parameters parameters;

    Vector<Vector3D> direction;
    Vector<Real>     weight;
    Vector<Size>     antipod;

    Matrix<Vector3D> m_direction;
    Matrix<Real>     m_weight;

    void read  (const Io& io);
    void write (const Io& io) const;

    accel inline Vector3D get_direction (
        const Size o,
        const Size r ) const;

    accel inline Real get_weight (
        const Size o,
        const Size r ) const;

    void set_1D_adaptive_rays (
        const Vector<Vector3D>& position);

    Vector<Vector3D> get_dirs (const Size o) const
    {
        cout << "Hey!" << endl;

        Vector<Vector3D> dirs;
                         dirs.resize(parameters.nrays());

        for (Size r = 0; r < parameters.nrays(); r++)
        {
            cout << "in" << endl;
            dirs[r] = m_direction(o,r);
            printf("%.16le,  %.16le\n", dirs[r].x(), dirs[r].y());
        }

        cout << "Out!" << endl;

        return dirs;
    }
};


#include "rays.tpp"
