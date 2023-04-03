#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/thermodynamics/temperature/temperature.hpp"


struct Frequencies
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    Matrix<Real> nu;                          ///< [Hz] frequencies (ordered in f) (p,f)
    Vector<Size> corresponding_line;          ///< [Hz] corresponding line to each frequency (f)

    Bool1 appears_in_line_integral;           ///< True if the frequency appears in line integral
    Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
    Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
    Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency

    //For imaging
    Size n_image_freqs;//set when imaging
    //Everything above will be resized when computing images (using Frequencies::resize_data)


    Frequencies (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;
    void resize_data(Size Nfreqs);



//    Size nbins = 0;    ///< number of extra bins per line
//    Size ncont = 0;    ///< number of background bins
};
