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
    vector <Size> nb_points_at_lvl;//keeps the number of points at each coarsening level
    Size current_nb_points;//is the number of points at the coarsest level

    vector<vector<bool>> mask_list;//stores the list of masks
    //true means point is still in grid, false means point is deleted

    //debug stuff for point deletion
    vector<std::map<Size, std::set<Size>>> reduced_neighbors_before;//without curr_point
    vector<vector<vector<Size>>> added_lines;//the lines in the order in which they are added
    vector<vector<vector<vector<Size>>>> added_tetras;//the tetrahedra in the order in which they are added (probably does not include all possible tetrahedra); might contain duplicates
    vector<Size> deleted_points;
    vector<std::map<Size, std::set<Size>>> reduced_neighbors_after;//without curr_point obviously
    bool debug_mode=false;



    void read  (const Io& io);
    void write (const Io& io) const;

    inline double orientation(const vector<Size> &plane, Size point);
    inline double calc_power(const vector<Size> &triangle, Size point);
    inline double calc_diff_abundance_with_neighbours(Size point, Size next_coars_lvl);
    inline void generate_initial_ears(const vector<Size> &neighbors_of_point, const vector<Size> &plane, std::map<Size, std::set<Size>> &neighbor_map,
     std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point);
    //inline void generate_new_ears(const std::set<vector<Size>> &neighbor_lines, vector<Size> &new_line,std::map<Size,std::set<Size>> &neighbor_map, std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point);
    inline void generate_new_ears(const vector<Size> &neighbors_of_point, const vector<Size> &plane, std::map<Size, std::set<Size>> &neighbor_map,
     std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point, Size &orient_point, std::set<std::vector<Size>> &forbidden_planes, std::set<std::vector<Size>> &edge_planes);

    inline void coarsen_grid(float perc_points_deleted);

    inline std::set<vector<Size>> calc_all_tetra_with_point(Size point, Size coars_lvl);
    inline Eigen::Vector<double,4> calc_barycentric_coords(const vector<Size> &triangle, Size point);
    //TODO:inline void rerefine_grid();
    inline vector<double> interpolate_vector(Size coarser_lvl, Size finer_lvl, const vector<double> &to_interpolate);
    inline void reset_grid();


    void read  ()       {read  (IoPython ("hdf5", parameters.model_name()));};
    void write () const {write (IoPython ("hdf5", parameters.model_name()));};

    int compute_inverse_line_widths     ();
    int compute_spectral_discretisation ();
    int compute_spectral_discretisation (const Real width);
    int compute_LTE_level_populations   ();
    int compute_radiation_field         ();
    int compute_Jeff                    ();
    int compute_level_populations       (
        // const Io   &io,
        const bool  use_Ng_acceleration,
        const long  max_niterations     );

    Double1 error_max;
    Double1 error_mean;
};

#include "model.tpp"
