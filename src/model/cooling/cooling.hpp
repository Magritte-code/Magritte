#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


// will contain structures with cooling rates
// should probably also include all/some methods to compute the cooling rates. (silly method of just taking the difference in intensity might not be too viable)

struct Cooling
{
    // // No longer a subclass of Model, so no direct access to the model parameters
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model
    Vector<Real> cooling_rate;//< for each point, contains the cooling rate
    //TODO: add individual cooling rate due to species?

    /// err, this is a data class; so no actual computation needs to take place here...
    // Vector<Real> mean_chiI;//< contains the angle_integrate mean intensity times opacity for each point, necessary for the flux version of the cooling computation
    // //TODO: think about removing mean_chiI, as the cooling rate can as well be derived using ∑_l χ_l J_l, in which χ_l is the integrated line opacity and J_l is the line integrated opacity
    // // in that case, one can easily separate the cooling introduced by each line
    // //TODO: think about defining the cooling rate per line, in addition to total cooling rate
    // Vector<Real> mean_grad_I;//< contains the (integrated over angle, frequency) intensity gradient for each point, necessary for the
    // Matrix<Real> line_grad_I;//< contains the (integrated over angle, frequency) line intensity gradient (=intensity gradient times line profile function, divided by opacity)
    // //After the feautrier refactor, add options to choose which things need to be computed

    Cooling (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    // void setup (Model& model);
    void write (const Io& io) const;

    // inline void compute_cooling_collisional(Model& model);
    // inline void compute_cooling_flux(Model& model);
    // inline void compute_cooling_grad_I(Model& model);

};
