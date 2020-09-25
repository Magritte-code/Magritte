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
    Double1 density_diff_at_point;//IMPORTANT find better name for this
    vector <Neighbors> neighbors_lists;
    Size current_nb_points;


    void read  (const Io& io);
    void write (const Io& io) const;


    inline double calc_diff_abundance_with_neighbours(int point, int next_coars_lvl);
    inline void coursen_grid(const float perc_points_deleted=0.5);
        //TODO
    inline void rerefine_grid();
        //TODO
    inline void reset_grid();


};


#include "model.tpp"
