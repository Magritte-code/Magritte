#pragma once


/////////////////////////////////////////
/// To Do: avoid having to define this...
/////////////////////////////////////////
#define MAX_NTHREADS 99
/////////////////////////////////////////

#define _USE_MATH_DEFINES
#include <cmath>

#include "../configure.hpp"
const string magritte_folder = MAGRITTE_FOLDER;

// Numerical constants
const double PI                         =     M_PI;         // pi
const double FOUR_PI                    = 4.0*M_PI;         // 4*pi
const double INVERSE_SQRT_PI            = 0.5*M_2_SQRTPI;   // 1/sqrt(pi)

const double ONE_THIRD                  = 1.0/3.0;
const double TWO_THIRDS                 = 2.0/3.0;
const double ONE_SIXTH                  = 1.0/6.0;

// Physical constants
const double CC                         = 2.99792458E+8;    ///< [m/s] speed of light
const double HH                         = 6.62607004E-34;   ///< [J*s] Planck's constant
const double KB                         = 1.38064852E-23;   ///< [J/K] Boltzmann's constant
const double AMU                        = 1.66053904E-27;   ///< [kg] proton mass
const double T_CMB                      = 2.7254800;        ///< [K] CMB temperature

const double SECONDS_IN_YEAR            = 3.1556926E+7;    // number of seconds in one year

const double CC_SQUARED                 = 8.98755179E+16;   // [m^2/s^2] speed of light squared

const double TWO_KB_OVER_AMU_CC_SQUARED = 2.0 * KB / (AMU * CC_SQUARED);   // 2.0*Kb/(AMU*c^2)
const double     HH_OVER_CC_SQUARED     = HH / CC_SQUARED;
const double TWO_HH_OVER_CC_SQUARED     = 2.0 * HH_OVER_CC_SQUARED;

const double HH_OVER_KB                 = HH / KB;
const double HH_OVER_FOUR_PI            = HH / FOUR_PI;
