#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

NEW_SYSTTOOLS_EXCEPT(invalid_SPPTpiReweight_Range);

namespace nusyst {

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

inline double GetSPPTpiReweight(double Tpi_GeV, double parameter_value){

  static double const P0 = 1.337359;
  static double const P0Err = 0.096961;
  static double const P1 = -2.769901;
  static double const P1Err = 0.318157;

  double X = Tpi_GeV;
  if(X>0.350) X = 0.350;

  double this_rw = P0 + P1 * X;
  if(this_rw<0) this_rw = 1.;

  //std::cout << "[GetSPPTpiReweight] Tpi_GeV = " << Tpi_GeV << ", RW = " << this_rw << std::endl;

  // dial = 0 : 1
  // dial = 1 : this_rw

  return (1-parameter_value) * 1. + parameter_value * this_rw;

}

inline double GetSPPTpiReweightMINERvA(double Tpi_GeV, double parameter_value){

  double X = Tpi_GeV*1000.; // GeV to MeV

  double this_rw = 1.;
  if( X < 10.000000) this_rw = 0.267183;
  else if( X >= 10.000000 && X < 15.000000) this_rw = 0.218322;
  else if( X >= 15.000000 && X < 20.000000) this_rw = 0.372796;
  else if( X >= 20.000000 && X < 25.000000) this_rw = 0.587210;
  else if( X >= 25.000000 && X < 30.000000) this_rw = 0.767524;
  else if( X >= 30.000000 && X < 36.000000) this_rw = 0.880305;
  else if( X >= 36.000000 && X < 42.000000) this_rw = 0.669767;
  else if( X >= 42.000000 && X < 48.000000) this_rw = 0.817111;
  else if( X >= 48.000000 && X < 54.000000) this_rw = 1.092730;
  else if( X >= 54.000000 && X < 60.000000) this_rw = 0.995627;
  else if( X >= 60.000000 && X < 66.000000) this_rw = 0.916708;
  else if( X >= 66.000000 && X < 72.000000) this_rw = 1.243540;
  else if( X >= 72.000000 && X < 78.000000) this_rw = 1.211460;
  else if( X >= 78.000000 && X < 84.000000) this_rw = 1.121870;
  else if( X >= 84.000000 && X < 90.000000) this_rw = 1.253250;
  else if( X >= 90.000000 && X < 96.000000) this_rw = 1.191510;
  else if( X >= 96.000000 && X < 102.000000) this_rw = 1.038230;
  else if( X >= 102.000000 && X < 110.000000) this_rw = 1.237920;
  else if( X >= 110.000000 && X < 125.000000) this_rw = 1.190560;
  else if( X >= 125.000000 && X < 140.000000) this_rw = 1.229080;
  else if( X >= 140.000000 && X < 155.000000) this_rw = 0.988201;
  else if( X >= 155.000000 && X < 175.000000) this_rw = 1.032940;
  else if( X >= 175.000000 && X < 200.000000) this_rw = 0.901374;
  else if( X >= 200.000000 && X < 225.000000) this_rw = 0.757748;
  else if( X >= 225.000000 && X < 250.000000) this_rw = 0.755932;
  else if( X >= 250.000000 && X < 275.000000) this_rw = 0.638574;
  else if( X >= 275.000000 && X < 300.000000) this_rw = 0.493987;
  else if( X >= 300.000000 && X < 325.000000) this_rw = 0.391947;
  else if( X >= 325.000000 && X < 350.000000) this_rw = 0.323265;
  else if( X >= 350.000000 && X < 400.000000) this_rw = 0.452765;
  else if( X >= 400.000000 && X < 500.000000) this_rw = 0.594541;
  else if( X >= 500.000000 && X < 700.000000) this_rw = 0.768459;
  else if( X >= 700.000000 && X < 1000.000000) this_rw = 0.658024;
  else this_rw = 0.873622;

  return (1-parameter_value) * 1. + parameter_value * this_rw;

}


};

