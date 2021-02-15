#pragma once

#include "model/mgController/mgController.hpp"
#include "model/mgController/naiveMG/naiveMG.hpp"


//The Multigrid controller helper structure
//just for the practical implementation
// for getting where we currently are in the multigrid sequence
struct MgControllerHelper : virtual public MgController
{
  private:
    enum class WhichImplementation {
      None,//Default type; no implementation
      NaiveMG
    };

    WhichImplementation current_implementation=WhichImplementation::None;

    //references to the different implementations
    NaiveMG naiveMGref;
    //TODO add many more

  public:

    //Inherited from mgController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   finish,                 //the entire multigrid procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    // };

    MgControllerHelper(){current_implementation=WhichImplementation::None;}

    //initializes the mgController
    //MgController(Size nb_levels, Size nb_pre_interpolation_steps);//TODO add much more

    //Because we cannot change constructor names, we use the default constructor just need to 'construct' using some regular functions
    inline void UseNaiveMG(Size nb_levels, Size finest_lvl);

    //returns the next action and updates what to do next
    Actions get_next_action();

    //returns the current level
    Size get_current_level();

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    void converged_on_current_grid();

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
};

#include "mgControllerHelper.tpp"
