#pragma once

#include "model/mrController/mrController.hpp"
//The multiresolution controller structure
// for getting where we currently are in the multiresolution sequence
struct NaiveMG : virtual public MrController
{
    private:

    Size max_n_iterations;          //the maximum number of iterations per level
    Size current_n_iterations = 0;  //the number of iterations done on the current grid

    Size finest_lvl           = 0;  //denotes the finest level the multiresolution sequence should go to
    Size max_level            = 0;  //number of levels going from 0 (finest grid) to n_levels-1 (coarsest grid)
    Size current_level        = 0;  //the current level

    bool is_next_action_set;        //checks whether the next action has been set
    Actions next_action;            //the next action (if it has been set)


    public:

    // Inherited from mrController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   stay_and_interpolate,   //remain on the same grid level, but also interpolate if necessary
    //   finish,                 //the entire multiresolution procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    //     }

    //Default necessary for mrControllerHelper
    // inline NaiveMG()=default;

    //Initializes the NaiveMG mrController
    inline NaiveMG(Size n_levels, Size finest_lvl, Size max_n_iterations);

    //returns the next action and updates what to do next
    inline Actions get_next_action() override;

    //returns the current level
    inline Size get_current_level() override;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    inline void converged_on_current_grid() override;

};


#include "naiveMG.tpp"
