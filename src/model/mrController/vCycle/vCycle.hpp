#pragma once

#include "model/mrController/mrController.hpp"
//The multiresolution controller structure
// for getting where we currently are in the multiresolution sequence
struct VCycle : virtual public MrController
{
    private:
    bool going_coarser;             //true when currently going to coarser grids
    bool first_upward     = true;   //first time going up should use standard interpolation of the level populations
    bool not_yet_iterated = true;   //checks whether we have already done a single iteration on the grid

    Size finest_lvl                 = 0;  //denotes the finest level the multiresolution sequence should go to
    Size max_level                  = 0;  //number of levels going from 0 (finest grid) to n_levels-1 (coarsest grid)
    Size current_level              = 0;  //the current level
    Size n_pre_interpolation_steps  = 1;  //number of iterations on the coarsest grid before interpolating

    Size max_n_iterations;          //the maximum number of cycles allowed for this multiresolution scheme
    Size current_n_iterations = 0;  //the current number of full cycles done

    bool is_next_action_set = false;  //checks whether the next action has been set
    Actions next_action;              //the next action (if it has been set)


    public:

    // Inherited from mrController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   finish,                 //the entire multiresolution procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    //     }

    //Default necessary for mrControllerHelper
    // inline VCycle()=default;

    //initializes the mrController
    inline VCycle(Size n_levels, Size finest_lvl, Size n_pre_interpolation_steps, Size max_n_iterations);//TODO add much more

    //returns the next action and updates what to do next
    inline Actions get_next_action() override;

    //returns the current level
    inline Size get_current_level() override;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    inline void converged_on_current_grid() override;

};


#include "vCycle.tpp"
