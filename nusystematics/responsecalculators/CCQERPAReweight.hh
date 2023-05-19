#ifndef nusystematics_RESPONSE_CALCULATORS_CCQERPAReweight_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCQERPAReweight_HH_SEEN

#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/exceptions.hh"

#include <ostream>

NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_RPA_tweak);

namespace nusyst {
  
class CCQERPAReweight
    : private nusyst::TemplateResponseCalculatorBase<3, false> {

  enum bin_indices { kIndex_enu = 0, kIndex_q0 = 1, kIndex_q3 = 2 };

public:

  CCQERPAReweight(fhicl::ParameterSet const &InputManifest) {
    LoadInputHistograms(InputManifest);
  }

  virtual bin_it_t GetBin(std::array<double, 3> const &kinematics) const {

    TH3 *firstHist = BinnedResponses.begin()->second.get();

    Int_t XBin = firstHist->GetXaxis()->FindFixBin(kinematics[kIndex_enu]);
    Int_t YBin = firstHist->GetYaxis()->FindFixBin(kinematics[kIndex_q0]);
    Int_t ZBin = firstHist->GetZaxis()->FindFixBin(kinematics[kIndex_q3]);
    // Hold events outside of the Valencia calculation phase space at the
    // closest valid bin.
    if (IsFlowBin(firstHist->GetXaxis(), XBin)) {
      XBin = (XBin == 0) ? XBin + 1 : XBin - 1;
    }
    if (IsFlowBin(firstHist->GetYaxis(), YBin)) {
      YBin = (YBin == 0) ? YBin + 1 : YBin - 1;
    }
    if (IsFlowBin(firstHist->GetZaxis(), ZBin)) {
      ZBin = (ZBin == 0) ? ZBin + 1 : ZBin - 1;
    }

#ifdef CCQERPAReweight_DEBUG
    printf("[CCQERPAReweight] (XBin, YBin, ZBin) = (%d, %d, %d)\n",XBin, YBin, ZBin);

#endif

    return firstHist->GetBin(XBin, YBin, ZBin);
  }

  double GetRPAReweight(double enu_GeV, double q0_GeV, double q3_GeV){

    double weight = 1;

    int bin3d = GetBin(std::array<double, 3>{{enu_GeV, q0_GeV, q3_GeV}});
    weight = GetVariation(e2i(tweak), bin2d);

    return weight;
  }

  std::string GetCalculatorName() const { return "CCQERPAReweight"; }
};
} // namespace nusyst

#endif
