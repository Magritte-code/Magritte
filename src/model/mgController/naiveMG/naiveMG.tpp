#pragma once

// Initializes the multigrid controller
// Finest level is for debug purposes, restrict the finest level until we iterate
inline NaiveMG::NaiveMG(Size nb_levels, Size finest_lvl, Size max_nb_iterations)
{
  //FIXME: check if nb_levels>=2 // otherwise we only have a single grid
  //FIXME: check if nb_pre_interpolation_steps>=1
  max_level=nb_levels-1;
  current_level=nb_levels-1;
  std::cout<<"current level:"<<current_level<<std::endl;
  this->finest_lvl=finest_lvl;

  //first, we go to the coarsest grid
  is_next_action_set=true;
  next_action=Actions::goto_coarsest;
  this->max_nb_iterations=max_nb_iterations;

}


inline Size NaiveMG::get_current_level()
{
  return current_level;
}

///Returns the next action and updates what to do next
inline MgController::Actions NaiveMG::get_next_action()
{
  //If converged, then use whatever recommended by converged_on_current_grid
  if (is_next_action_set)
  {
    is_next_action_set=false;
    std::cout<<"Returning next action set"<<std::endl;
    return next_action;
  }else{//otherwise just remain on the current grid, iterating until convergence/max_nb_iterations reached
    if (current_nb_iterations<max_nb_iterations)
    {
      std::cout<<"Returning stay"<<std::endl;
      current_nb_iterations++;
      return Actions::stay;
    }else{
      if (current_level>finest_lvl)
      {
        std::cout<<"Returning interpolate_levelpops"<<std::endl;
        current_level--;
        current_nb_iterations=0;
        return Actions::interpolate_levelpops;
      }else{
        std::cout<<"Returning finish"<<std::endl;
        return Actions::finish;
      }
    }
  }
}

//Call this when the solution is converged on the current grid.
//Sets the state such that the next action will be something else than stay (skips the following 'stay's)
inline void NaiveMG::converged_on_current_grid()
{
  if (current_level==finest_lvl)
  {
    next_action=Actions::finish;
    is_next_action_set=true;
    return;
  }
  else if (current_level>finest_lvl)
  {
    //reset current_nb_iterations, as we interpolate to the finer grid
    current_nb_iterations=0;
    next_action=Actions::interpolate_levelpops;
    is_next_action_set=true;
    std::cout<<"current level: "<<current_level<<std::endl;
    current_level--;
    std::cout<<"finest level: "<<finest_lvl<<std::endl;
    return;
  }
  else
  {
    next_action=Actions::finish;
    //TODO also return error
    std::cout<<"Somehow, you are currently on a level finer than the finest level you allowed"<<std::endl;
    std::cout<<"Finishing the multigrid computation either way"<<std::endl;
    std::cout<<"Finest level: "<<finest_lvl<<"Current level: "<<current_level<<std::endl;
    return;
  }
}
