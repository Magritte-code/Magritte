#pragma once


#include "io/io.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"
#include "radiation/radiation.hpp"


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    // Thermodynamics thermodynamics;
    // Lines          lines;
    // Radiation      radiation;


    int curr_coarsening_lvl=0;
    int max_reached_coarsening_lvl=0;
    std::multimap<int,double> density_diff_map;//stores (point,density_diff)
    std::multimap<double,int> rev_density_diff_map;//for ease of accessing data, i also store the reverse map
      //@Frederik: if you know a better data type to store this, let me know
    vector <Neighbors> neighbors_lists;
    Size current_nb_points;


    void read  (const Io& io);
    void write (const Io& io) const;

    inline double calc_power(vector<int> triangle, int point);
    inline double calc_diff_abundance_with_neighbours(int point, int next_coars_lvl);
    inline void coarsen_grid(const float perc_points_deleted=0.5);
        //TODO
    inline void rerefine_grid();
    inline void reset_grid();


};


#include "model.tpp"
