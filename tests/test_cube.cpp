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

// Vector<Vector3D> positions= Vector<Vector3D>(9);
// positions.resize(9);
// positions[0] = Vector3D(0,1,1);
// positions[1] = Vector3D(1,1,1);
// positions[2] = Vector3D(0,0,1);
// positions[3] = Vector3D(1,0,1);
// positions[4] = Vector3D(0,1,0);
// positions[5] = Vector3D(1,1,0);
// positions[6] = Vector3D(0,0,0);
// positions[7] = Vector3D(1,0,0);
// positions[8] = Vector3D(1/2,1/2,1/2);


//TODO: also add coordinates...
// and boundary
// and ofcourse parameters.npoints


Model model;
model.parameters.set_npoints(9);
//vector<Size> temp_point2boundary{0,1,2,3,4,5,6,7,9};
model.geometry.boundary.point2boundary.resize(9);
//I still cant figure out how shallow copies work with paracabs...
//model.geometry.boundary.point2boundary=Vector<Size>(temp_point2boundary);
model.geometry.boundary.point2boundary[0]=0;
model.geometry.boundary.point2boundary[1]=1;
model.geometry.boundary.point2boundary[2]=2;
model.geometry.boundary.point2boundary[3]=3;
model.geometry.boundary.point2boundary[4]=4;
model.geometry.boundary.point2boundary[5]=5;
model.geometry.boundary.point2boundary[6]=6;
model.geometry.boundary.point2boundary[7]=7;
model.geometry.boundary.point2boundary[8]=9;


model.geometry.points.position.resize(9);
model.geometry.points.position[0]=Vector3D(0,1,1);
model.geometry.points.position[1]=Vector3D(1,1,1);
model.geometry.points.position[2]=Vector3D(0,0,1);
model.geometry.points.position[3]=Vector3D(1,0,1);
model.geometry.points.position[4]=Vector3D(0,1,0);
model.geometry.points.position[5]=Vector3D(1,1,0);
model.geometry.points.position[6]=Vector3D(0,0,0);
model.geometry.points.position[7]=Vector3D(1,0,0);
model.geometry.points.position[8]=Vector3D(1./2,1./2,1./2);

std::cout<<"Position of middle point: "<<model.geometry.points.position[8].x() << ","<<model.geometry.points.position[8].y() << ","<<model.geometry.points.position[8].z() << std::endl;

model.geometry.points.curr_neighbors.set_all_neighbors(n_neighbors, neighbors);



model.chemistry.species.abundance=Double2(9,std::vector<double>(2,0));

std::cout << model.parameters.get_npoints() << std::endl;

for (Size i = 0; i < model.parameters.get_npoints(); ++i)
    std::cout << model.geometry.boundary.point2boundary[i] << ' ';
std::cout << std::endl;
//should print 0,1,2,3,4,5,6,7,9

model.coarsen_grid(0.2);//if i am correct, this should round down on the nb of points deleted
for (auto i = model.neighbors_lists[0].n_neighbors.begin(); i != model.neighbors_lists[0].n_neighbors.end(); ++i)
    std::cout << *i << ' ';
std::cout << std::endl;
//printing out a vector
for (auto i = model.neighbors_lists[1].n_neighbors.begin(); i != model.neighbors_lists[1].n_neighbors.end(); ++i)
    std::cout << *i << ' ';
std::cout << std::endl;

for (Size i = 0; i<model.parameters.get_npoints(); ++i){
  for(Size j = 0; j<model.neighbors_lists[0].n_neighbors[i]; ++j){
    std::cout << model.neighbors_lists[0].neighbors[i][j];}
    std::cout << ',';}
std::cout << std::endl;

for (Size i = 0; i<model.parameters.get_npoints(); ++i){
  for(Size j = 0; j<model.neighbors_lists[1].n_neighbors[i]; ++j){
    std::cout << model.neighbors_lists[1].neighbors[i][j];}
    std::cout << ',';}
std::cout << std::endl;

}
