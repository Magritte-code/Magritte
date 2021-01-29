// Voronoi calculation example code
//
// Author   : Chris H. Rycroft (Harvard University / LBL)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

// adapted for testing purposes

#include "voro++.hh"
#include <iostream>
#include <fstream>
using namespace voro;

// Set up constants for the container geometry
const double x_min=-2.5,x_max=2.5;
const double y_min=-2.5,y_max=2.5;
const double z_min=-1.1,z_max=1.1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
const int n_x=10,n_y=10,n_z=10;

// Set the number of particles that are going to be randomly introduced
const int particles=9;
// const int particles=200;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	double x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

  // std::vector<double> xvals{0,-1,1,0,-1,1,0,-1,1,2,-2,0,0};
  // std::vector<double> yvals{0,0,0,1,1,1,-1,-1,-1,0,0,2,-2};
  // std::vector<double> zvals{0,0,0,0,0,0,0,0,0,0,0,0,0};
	std::vector<double> xvals{0,-1,1,-1,1,2,-2,0,0};
	std::vector<double> yvals{0,1,1,-1,-1,0,0,2,-2};
	std::vector<double> zvals{0,0,0,0,0,0,0,0,0};

  particle_order p_order;

  voronoicell_neighbor cell;

	// add particles into the container
	for(i=0;i<particles;i++) {
		// x=x_min+rnd()*(x_max-x_min);
		// y=y_min+rnd()*(y_max-y_min);
		// z=z_min+rnd()*(z_max-z_min);
		con.put(p_order,10+i,xvals[i],yvals[i],zvals[i]);
	}

  c_loop_order l_order(con,p_order);

  //back to the beginning
  l_order.start();
  int idx=0;
  //use inc() for incrementing (returns false when not finding a next one)
  std::cout<<"p_order size: "<<p_order.size<<std::endl;
  do {
  con.compute_cell(cell,l_order);
  std::vector<int> test_neighbors;
  cell.neighbors(test_neighbors);
  std::cout<<"neighbor of cell "<<l_order.pid()<<": ";
  //note:negative values refer to the boundaries
  for (int i:test_neighbors)
  {
    std::cout<<i<<", ";
  }
  std::cout<<std::endl;
  idx++;
  }
  while(l_order.inc());

	// Sum up the volumes, and check that this matches the container volume
	// double vvol=con.sum_cell_volumes();
	// printf("Container volume : %g\n"
	//        "Voronoi volume   : %g\n"
	//        "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");
  // use gnuplot with: splot "random_points_v.gnu" with lines
  // in order to see something
	const char *vars = "%n";
	con.print_custom(vars,"test.txt");

}
