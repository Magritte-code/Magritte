#pragma once

#include "model/mgController/mgController.hpp"
//The Multigrid controller structure
// for getting where we currently are in the multigrid sequence
struct VCycle : virtual public MgController
{
  private:
    bool going_coarser;// true when currently going to coarser grids
    bool first_upward=true;//first time going up should use standard interpolation of the level populations
    bool not_yet_iterated=true;//checks whether we have already done a single iteration on the grid

    // Size nb_iterations_on_current_grid_remaining;// number of iterations remaining on current grid
    Size finest_lvl=0;
    Size max_level=0;//number of levels going from 0 (finest grid) to nb_levels-1 (coarsest grid)
    Size current_level=0;//the current level
    //Size min_level_visited=0;
    Size nb_pre_interpolation_steps=1;//number of iterations on the current grid before interpolating/coarsening

    Size max_nb_iterations;
    Size current_nb_iterations=0;

    bool is_next_action_set=false;//checks whether the next action has been set
    Actions next_action;//the next action (if it has been set)
    //for now, I implement the V-cycle


  public:

    // Inherited from mgController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   finish,                 //the entire multigrid procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    //     }

    //Default necessary for mgControllerHelper
    // inline VCycle()=default;

    //initializes the mgController
    inline VCycle(Size nb_levels, Size finest_lvl, Size nb_pre_interpolation_steps, Size max_nb_iterations);//TODO add much more

    //returns the next action and updates what to do next
    inline Actions get_next_action() override;

    //returns the current level
    inline Size get_current_level() override;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    inline void converged_on_current_grid() override;

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
};


#include "vCycle.tpp"
