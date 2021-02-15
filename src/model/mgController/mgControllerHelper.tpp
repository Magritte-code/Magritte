#pragma once

// #include "model/mgController/naiveMG/naiveMG.hpp"

    //Because we cannot change constructor names, we use the default constructor just need to 'construct' using some regular functions
inline void MgControllerHelper::UseNaiveMG(Size nb_levels, Size finest_lvl)
{
  naiveMGref=NaiveMG(nb_levels, finest_lvl);
  current_implementation=WhichImplementation::NaiveMG;
}

    //returns the next action and updates what to do next
MgController::Actions MgControllerHelper::get_next_action()
{
  switch (current_implementation) {
    case WhichImplementation::NaiveMG:
      return naiveMGref.get_next_action();
      break;
    case WhichImplementation::None:
    default:
    //TODO: raise error
      break;
  }
  return MgController::Actions::error;
}

//returns the current level
Size MgControllerHelper::get_current_level()
{
  switch (current_implementation) {
    case WhichImplementation::NaiveMG:
      return naiveMGref.get_current_level();
      break;
    case WhichImplementation::None:
    default:
    //TODO: raise error
      break;
  }
  return 0;
}

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
void MgControllerHelper::converged_on_current_grid()
{
  switch (current_implementation) {
    case WhichImplementation::NaiveMG:
      naiveMGref.converged_on_current_grid();
      break;
    case WhichImplementation::None:
    default:
    //TODO: raise error
      break;
  }
}

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
