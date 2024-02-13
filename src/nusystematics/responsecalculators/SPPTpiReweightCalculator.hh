#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

NEW_SYSTTOOLS_EXCEPT(invalid_SPPTpiReweight_Range);

namespace nusyst {

inline double GetSPPTpiReweight(double Tpi_GeV, double parameter_value){


  static double const P0 = 1.347965;
  static double const P0Err = 0.097783;
  static double const P1 = -2.817372;
  static double const P1Err = 0.320814;

  double X = Tpi_GeV;
  if(X>0.350) X = 0.350;

  double this_rw = P0 + P1 * X;
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
  if( X < 0.025000) this_rw = 1.253255;
  else if( X >= 0.025000 && X < 0.050000) this_rw = 1.589738;
  else if( X >= 0.050000 && X < 0.100000) this_rw = 1.733869;
  else if( X >= 0.100000 && X < 0.200000) this_rw = 1.651728;
  else if( X >= 0.200000 && X < 0.300000) this_rw = 1.659705;
  else if( X >= 0.300000 && X < 0.400000) this_rw = 1.584229;
  else if( X >= 0.400000 && X < 0.500000) this_rw = 1.703793;
  else if( X >= 0.500000 && X < 0.700000) this_rw = 1.475510;
  else if( X >= 0.700000 && X < 1.000000) this_rw = 1.456727;
  else if( X >= 1.000000 && X < 1.300000) this_rw = 1.252215;
  else if( X >= 1.300000 && X < 2.000000) this_rw = 1.048199;
  else if( X >= 2.000000 && X < 3.000000) this_rw = 1.650489;
  else{
    this_rw = 1.;
  }

  return (1-parameter_value) * 1. + parameter_value * this_rw;

}

};

