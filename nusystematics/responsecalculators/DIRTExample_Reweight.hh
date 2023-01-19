#ifndef nusystematics_RESPONSE_CALCULATORS_DIRTEXAMPLE_REWEIGHT_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_DIRTEXAMPLE_REWEIGHT_HH_SEEN

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

  inline double GetXSecCV(double q3, double q0){
    if(q3>0.5) return 1.;
    if(q0<0) return 1.; // useless line to avoid unused-parameter error
    return 1.;
  }
  inline double GetXSecAltModel(double q3, double q0){
    if(q3>0.5) return 0.;
    if(q0<0) return 0.; // useless line to avoid unused-parameter error
    else return 1.1;
  }

  inline double GetDIRTExampleWeight(double q3, double q0, 
                                     double parameter_value) {

    // parameter_value = 0 : CV
    // parameter_value = 1 : AltModel
    double xsec_cv = GetXSecCV(q3, q0);
    double xsec_alt = GetXSecAltModel(q3, q0);
    double new_xsec = parameter_value*xsec_alt + (1.-parameter_value)*xsec_cv;
    return new_xsec/xsec_cv;

  }

} // namespace nusyst

#endif
