#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "model/thermodynamics/temperature/temperature.hpp"
#include "tools/types.hpp"

struct Frequencies {
    std::shared_ptr<Parameters> parameters; ///< data structure containing model parameters

    Matrix<Real> nu;                 ///< [Hz] frequencies (ordered in block of lines) (p,f)
    Vector<Size> corresponding_line; ///< corresponding line to each frequency (f)
    Matrix<Real> sorted_nu;          ///< [Hz] sorted frequencies (ordered in f at each point) (p,f)
    Matrix<Size> corresponding_nu_index; ///< corresponding (original) frequency index for each
                                         ///< frequency at each point (p,f)

    Bool1 appears_in_line_integral; ///< True if the frequency appears in line integral
    Size1 corresponding_l_for_spec; ///< number of line species corresponding to frequency
    Size1 corresponding_k_for_tran; ///< number of transition corresponding to frequency
    Size1 corresponding_z_for_line; ///< number of line number corresponding to frequency

    Frequencies(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;
    void resize_data(Size Nfreqs);

    //    Size nbins = 0;    ///< number of extra bins per line
    //    Size ncont = 0;    ///< number of background bins
};
