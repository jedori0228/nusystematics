#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

NEW_SYSTTOOLS_EXCEPT(invalid_SPPTpiReweight_Range);

namespace nusyst {

inline double GetSPPTpiReweight(double Tpi_GeV, double parameter_value){

  static double NormSF = 2.17029;

  static double const P0 = 1.0000000;
  static double const P0Err = 0.;
  static double const P1 = -2.051462;
  static double const P1Err = 0.150963;

  double X = Tpi_GeV;
  if(X>0.350) X = 0.350;

  double this_rw = NormSF * (P0 + P1 * X);
  if(this_rw<0) this_rw = 1.;

  //std::cout << "[GetSPPTpiReweight] Tpi_GeV = " << Tpi_GeV << ", RW = " << this_rw << std::endl;

  // dial = 0 : 1
  // dial = 1 : this_rw

  return (1-parameter_value) * 1. + parameter_value * this_rw;

}

inline double GetSPPQ2Reweight(double Q2_GeV2, double parameter_value){

  double X = Q2_GeV2;
  if(X>=3.000000) X = 3.000000;

  double this_rw = 1.;
  if( Q2_GeV2 < 0.025000) this_rw = 1.143409;
  else if( X >= 0.025000 && Q2_GeV2 < 0.050000) this_rw = 1.384572;
  else if( X >= 0.050000 && Q2_GeV2 < 0.100000) this_rw = 1.523289;
  else if( X >= 0.100000 && Q2_GeV2 < 0.200000) this_rw = 1.504776;
  else if( X >= 0.200000 && Q2_GeV2 < 0.300000) this_rw = 1.530361;
  else if( X >= 0.300000 && Q2_GeV2 < 0.400000) this_rw = 1.477627;
  else if( X >= 0.400000 && Q2_GeV2 < 0.500000) this_rw = 1.559800;
  else if( X >= 0.500000 && Q2_GeV2 < 0.700000) this_rw = 1.389166;
  else if( X >= 0.700000 && Q2_GeV2 < 1.000000) this_rw = 1.373498;
  else if( X >= 1.000000 && Q2_GeV2 < 1.300000) this_rw = 1.215221;
  else if( X >= 1.300000 && Q2_GeV2 < 2.000000) this_rw = 1.044024;
  else if( X >= 2.000000 && Q2_GeV2 < 3.000000) this_rw = 1.574511;
  else{
    this_rw = 1.;
  }

  return (1-parameter_value) * 1. + parameter_value * this_rw;

}

};

