#pragma once


#include "io/io.hpp"
#include "io/python/io_python.hpp"
#include "parameters/parameters.hpp"
#include "tools/types.hpp"
#include "geometry/geometry.hpp"
#include "chemistry/chemistry.hpp"
#include "thermodynamics/thermodynamics.hpp"
#include "lines/lines.hpp"
#include "radiation/radiation.hpp"
#include <set>


struct Model
{
    Parameters     parameters;
    Geometry       geometry;
    Chemistry      chemistry;
    Thermodynamics thermodynamics;
    Lines          lines;
    Radiation      radiation;

    enum SpectralDiscretisation {None, SD_Lines, SD_Image}
         spectralDiscretisation = None;

    Model () {};
    Model (const string name)
    {
        parameters.model_name() = name;
        read ();
    }


    Size curr_coarsening_lvl=0;
    Size max_reached_coarsening_lvl=0;
    std::multimap<Size,double> density_diff_map;//stores (point,density_diff)
    std::multimap<double,Size> rev_density_diff_map;//for ease of accessing data, i also store the reverse map
      //@Frederik: if you know a better data type to store this, let me know
    vector <Neighbors> neighbors_lists;
    Size current_nb_points;

    vector<vector<bool>> mask_list;//stores the list of masks
    //true means point is still in grid, false means point is deleted


    void read  (const Io& io);
    void write (const Io& io) const;

    inline double calc_power(const vector<Size> &triangle, Size point);
    inline double calc_diff_abundance_with_neighbours(Size point, Size next_coars_lvl);
    inline void generate_new_ears(const vector<Size> &neighbors_of_point,const vector<Size> &plane,std::map<Size,std::set<Size>> &neighbor_map, std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point);
    inline void coarsen_grid(float perc_points_deleted);

    inline Eigen::Vector<double,4> calc_barycentric_coords(const vector<Size> &triangle, Size point);
    //TODO:inline void rerefine_grid();
    inline vector<double> interpolate_vector(Size coarser_lvl, Size finer_lvl)
    inline void reset_grid();


    void read  ()       {read  (IoPython ("hdf5", parameters.model_name()));};
    void write () const {write (IoPython ("hdf5", parameters.model_name()));};

    int compute_inverse_line_widths     ();
    int compute_spectral_discretisation ();
    int compute_spectral_discretisation (const Real width);
    int compute_LTE_level_populations   ();
    int compute_radiation_field         ();
    int compute_Jeff                    ();
};

#include "model.tpp"
