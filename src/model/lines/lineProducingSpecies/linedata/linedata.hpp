#pragma once


#include <Eigen/Core>

#include "io/io.hpp"
#include "tools/types.hpp"
#include "collisionPartner/collisionPartner.hpp"


class Linedata
{
    public:
        Size   num;                              ///< number of line producing species
        string sym;                              ///< symbol of line producing species
        Real   inverse_mass;                     ///< 1/mass of line producing species

        Size nlev;                               ///< number of levels
        Size nrad;                               ///< number of radiative transitions

        Size1 irad;                              ///< level index of radiative transition
        Size1 jrad;                              ///< level index of radiative transition

        Real1 energy;                            ///< energy of level
        Real1 weight;                            ///< weight of level (statistical)

        Real1 frequency;                         ///< frequency corresponding to each transition

        Real1 A;                                 ///< Einstein A  (spontaneous emission)
        Real1 Ba;                                ///< Einstsin Ba (absorption)
        Real1 Bs;                                ///< Einstein Bs (stimulated emission)

        Size ncolpar;                            ///< number of collision partners

        // Collision partners
        std::vector <CollisionPartner> colpar;   ///< Vector containing collision partner data

        Size ncol_tot;

        void read  (const Io &io, const Size l);
        void write (const Io &io, const Size l) const;
};
