#ifndef nusystematics_SYSTPROVIDERS_MINERvASPPReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_MINERvASPPReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/MINERvASPPReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class MINERvASPPReweight : public nusyst::IGENIESystProvider_tool {

public:

  NEW_SYSTTOOLS_EXCEPT(invalid_engine_state);

  explicit MINERvASPPReweight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~MINERvASPPReweight();

private:

  fhicl::ParameterSet tool_options;

  // Find and save pidx and use later when filling the reweights
  // 1) Q2
  // If you are using GENIEv3, G18 or AR23 or similar,
  // Q2 must be corrected prior to Tpi
  size_t pidx_Q2;
  // 2) Tpi correction
  size_t pidx_Tpi;

};

#endif
