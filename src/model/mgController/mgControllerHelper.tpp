#pragma once

#include <memory>

//just assign the properly initialize mgController to this.
inline MgControllerHelper::MgControllerHelper(std::shared_ptr<MgController> implementation_instance_ptr)
{
  // MgController* temp_mg_control=dynamic_cast<MgController*>(&implementation_instance);
  // implementation_instance_ptr=dynamic_cast<MgController*>(&implementation_instance);
  // implementation_instance_ptr=std::make_shared<MgController>(implementation_instance_ptr);
  this->implementation_instance_ptr=implementation_instance_ptr;
}

//Returns the next action and updates what to do next
MgController::Actions MgControllerHelper::get_next_action()
{ //check if not nullptr
  if (implementation_instance_ptr)
  {
  std::cout<<"Getting next action"<<std::endl;
  return (*implementation_instance_ptr).get_next_action();
  }
  else
  {
    std::cout<<"You forgot to properly initialize mgControllerHelper;"<<std::endl;
    std::cout<<"Did you forget to use setup_multigrid?"<<std::endl;
    throw std::runtime_error("Error: mgControllerHelper pointer not set");
    // return MgController::Actions::error;
  }
}

//Returns the current level // should be deprecated
Size MgControllerHelper::get_current_level()
{
  if (implementation_instance_ptr)
  {
    return (*implementation_instance_ptr).get_current_level();
  }
  else
  {
    std::cout<<"You forgot to properly initialize mgControllerHelper;"<<std::endl;
    std::cout<<"Did you forget to use setup_multigrid?"<<std::endl;
    throw std::runtime_error("Error: mgControllerHelper pointer not set");
    // return parameters.npoints() //if uncomment, do not forget to include parameters somewhere in the hpp file
  }
}

///Call this when the solution is converged on the current grid.
///Sets the state such that the next action will be something else than stay (skips the following 'stay's)
void MgControllerHelper::converged_on_current_grid()
{
  if (implementation_instance_ptr)
  {
    (*implementation_instance_ptr).converged_on_current_grid();
  }
  else
  {
    std::cout<<"You forgot to properly initialize mgControllerHelper;"<<std::endl;
    std::cout<<"Did you forget to use setup_multigrid?"<<std::endl;
    throw std::runtime_error("Error: mgControllerHelper pointer not set");
  }
}
