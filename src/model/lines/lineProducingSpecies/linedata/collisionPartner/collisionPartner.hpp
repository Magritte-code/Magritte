#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

struct CollisionPartner {
    Size num_col_partner;   ///< species number corresponding to collision partner
    string orth_or_para_H2; ///< stores whether it is ortho or para (if it is H2)

    Size ntmp; ///< number of defined temperatures
    Size ncol; ///< number of collisional transitions

    Size1 icol; ///< level index of collisional transition
    Size1 jcol; ///< level index of collisional transition

    Real1 tmp; ///< Collision temperatures for each partner

    Real2 Ce; ///< Collisional excitation rates for each temperature
    Real2 Cd; ///< Collisional de-excitation rates for each temperature

    Real1 Ce_intpld; ///< interpolated Collisional excitation
    Real1 Cd_intpld; ///< interpolated Collisional de-excitation

    void read(const Io& io, const Size l, const Size c);
    void write(const Io& io, const Size l, const Size c) const;

    inline void adjust_abundance_for_ortho_or_para(
        const Real temperature_gas, Real& abundance) const;

    inline void interpolate_collision_coefficients(const Real temperature_gas);
};

#include "collisionPartner.tpp"
