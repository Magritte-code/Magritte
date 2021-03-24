#include "boundary.hpp"


const string prefix = "geometry/boundary/";


void Boundary :: read (const Io& io)
{
    cout << "Reading boundary..." << endl;

    // Resize boundary
    point2boundary.resize (parameters.npoints());

    // Initialise
    for (Size p = 0; p < parameters.npoints(); p++)
    {
        point2boundary[p] = parameters.npoints();
    }

    // Read boundary list
    parameters.set_nboundary (io.get_length (prefix+"boundary2point"));

    boundary2point.resize (parameters.nboundary());

    io.read_list (prefix+"boundary2point", boundary2point);

    // Set helper variables to identify the boundary
    for (Size b = 0; b < parameters.nboundary(); b++)
    {
        point2boundary[boundary2point[b]] = b;
    }


    // Set boundary conditions
    boundary_condition  .resize (parameters.nboundary());
    boundary_temperature.resize (parameters.nboundary());

    Size1 boundary_condition_int (parameters.nboundary());

    io.read_list (prefix+"boundary_temperature", boundary_temperature);
    io.read_list (prefix+"boundary_condition", boundary_condition_int);


    for (Size b = 0; b < parameters.nboundary(); b++)
    {
        if (boundary_condition_int[b] == 0) boundary_condition[b] = Zero;
        if (boundary_condition_int[b] == 1) boundary_condition[b] = Thermal;
        if (boundary_condition_int[b] == 2) boundary_condition[b] = CMB;
    }


    boundary2point.copy_vec_to_ptr ();
    point2boundary.copy_vec_to_ptr ();

    boundary_condition  .copy_vec_to_ptr ();
    boundary_temperature.copy_vec_to_ptr ();
}


void Boundary :: write (const Io& io) const
{
    cout << "Writing boundary..." << endl;

    io.write_list (prefix+"boundary2point", boundary2point);

    Size1 boundary_condition_int (parameters.nboundary());

    for (Size b = 0; b < parameters.nboundary(); b++) switch (boundary_condition[b])
    {
        case Zero    : boundary_condition_int[b] = 0;
        case Thermal : boundary_condition_int[b] = 1;
        case CMB     : boundary_condition_int[b] = 2;
    }

    io.write_list (prefix+"boundary_temperature", boundary_temperature);
    io.write_list (prefix+"boundary_condition", boundary_condition_int);
}


BoundaryCondition Boundary :: set_boundary_condition (const Size b, const BoundaryCondition cd)
{
    boundary_condition.resize(parameters.nboundary());
    boundary_condition[b] = cd;

    return boundary_condition[b];
}


BoundaryCondition Boundary :: get_boundary_condition (const Size b) const
{
    return boundary_condition[b];
}
