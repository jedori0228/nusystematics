#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include "Physics/NuclearState/LocalFGM.h"

#include <cmath>

namespace nusyst {
    // Default WS params for Ar40
    const double kAr40Radius = 3.53;
    const double kAr40SkinDepth = 0.54;

    inline double GetRadiusFromKF(double KF, bool isProton) {
        // Calculate the radius of the nucleus from the local Fermi momentum
        //     KF = (3*pi^2 * n * rho)^(1/3) / hbarc
        //Assuming only Ar40 for now
        int A = 40;
        int N = 22;
        int Z = 18;
        
        double norm  = (3./(4.*genie::constants::kPi*TMath::Power(kAr40Radius,3)))*1./(1.+TMath::Power((genie::constants::kPi*kAr40SkinDepth/kAr40Radius),2));
    
        double n_nucleons = (isProton) ? Z : N;
        double log_term = (3*norm*n_nucleons*genie::constants::kPi2)/TMath::Power(KF*genie::units::fermi,3) - 1;
        if(log_term <= 0){
            return 0;
        }
        return kAr40Radius + kAr40SkinDepth*TMath::Log(log_term);
      }
    
      inline double GetWeightFomKF(double KF, double c_val, double z_val, bool isProton) {
        // Calculate the weight from the local Fermi momentum
        //     KF = (3*pi^2 * n * rho)^(1/3) / hbarc
        //     weight = 1 + parameter_value * (KF - KF0)
       
    
        double radius = GetRadiusFromKF(KF, isProton);
        radius = std::max(radius, 0.0); // No negative radii!

        double ref_prob_density = genie::utils::nuclear::Density(radius, 40);
        double new_prob_density = genie::utils::nuclear::DensityWoodsSaxon(radius, c_val, z_val);
    
        return new_prob_density / ref_prob_density;
    
      }
                                  

} // namespace nusyst


