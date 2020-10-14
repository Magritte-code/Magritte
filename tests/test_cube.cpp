#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"
#include "tools/types.hpp"

int main (int argc, char **argv)
{

vector<Size> n_neighbors{6,4,6,6,6,6,4,6,8};
//note: the last two are due to splitting the sides of the cube
vector<Size> neighbors {1,2,4,8, 3,5,
                        0,3,5,8,
                        0,3,6,8, 4,7,
                        1,2,7,8, 0,5,
                        0,5,6,8, 2,7,
                        1,4,7,8, 0,3,
                        2,4,7,8,
                        3,5,6,8, 2,4,
                        0,1,2,3,4,5,6,7};

Vector<Vector3D> positions= Vector<Vector3D>(9);
positions.resize(9);
positions[0] = Vector3D(0,1,1);
positions[1] = Vector3D(1,1,1);
positions[2] = Vector3D(0,0,1);
positions[3] = Vector3D(1,0,1);
positions[4] = Vector3D(0,1,0);
positions[5] = Vector3D(1,1,0);
positions[6] = Vector3D(0,0,0);
positions[7] = Vector3D(1,0,0);
positions[8] = Vector3D(1/2,1/2,1/2);


//TODO: also add coordinates...
// and boundary
// and ofcourse parameters.npoints


Model model;
model.geometry.points.curr_neighbors.set_all_neighbors(n_neighbors, neighbors);
model.parameters.set_npoints(9);


}
