#pragma once

#include <memory>
// #include "model/mgController/naiveMG/naiveMG.hpp"

    //Because we cannot change constructor names, we use the default constructor just need to 'construct' using some regular functions
// inline void MgControllerHelper::UseNaiveMG(Size nb_levels, Size finest_lvl)
// {
//   naiveMGref=NaiveMG(nb_levels, finest_lvl);
//   current_implementation=WhichImplementation::NaiveMG;
//
// }

//just assign the properly initialize mgController to this.
// template <typename MgImplementation>
inline MgControllerHelper::MgControllerHelper(std::shared_ptr<MgController> implementation_instance_ptr)
{
  // MgController* temp_mg_control=dynamic_cast<MgController*>(&implementation_instance);
  // implementation_instance_ptr=dynamic_cast<MgController*>(&implementation_instance);
  // implementation_instance_ptr=std::make_shared<MgController>(implementation_instance_ptr);
  this->implementation_instance_ptr=implementation_instance_ptr;
}

    //returns the next action and updates what to do next
// template <typename MgImplementation>
MgController::Actions MgControllerHelper::get_next_action()
{
  // switch (current_implementation) {
  //   case WhichImplementation::NaiveMG:
  //     return naiveMGref.get_next_action();
  //     break;
  //   case WhichImplementation::None:
  //   default:
  //   //TODO: raise error
  //     break;
  // }
  // return MgController::Actions::error;
  std::cout<<"Getting next action"<<std::endl;
  return (*implementation_instance_ptr).get_next_action();
}

//returns the current level
// template <typename MgImplementation>
Size MgControllerHelper::get_current_level()
{
  // switch (current_implementation) {
  //   case WhichImplementation::NaiveMG:
  //     return naiveMGref.get_current_level();
  //     break;
  //   case WhichImplementation::None:
  //   default:
  //   //TODO: raise error
  //     break;
  // }
  // return 0;
  return (*implementation_instance_ptr).get_current_level();
}

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
// template <typename MgImplementation>
void MgControllerHelper::converged_on_current_grid()
{
  // switch (current_implementation) {
  //   case WhichImplementation::NaiveMG:
  //     naiveMGref.converged_on_current_grid();
  //     break;
  //   case WhichImplementation::None:
  //   default:
  //   //TODO: raise error
  //     break;
  // }
  (*implementation_instance_ptr).converged_on_current_grid();
}

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
