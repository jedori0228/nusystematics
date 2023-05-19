#ifndef nusystematics_MAKE_INSTANCE_SEEN
#define nusystematics_MAKE_INSTANCE_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/systproviders/BeRPAWeight_tool.hh"
#include "nusystematics/systproviders/EbLepMomShift_tool.hh"
#include "nusystematics/systproviders/FSILikeEAvailSmearing_tool.hh"
#include "nusystematics/systproviders/GENIEReWeight_tool.hh"
#include "nusystematics/systproviders/MINERvAE2p2h_tool.hh"
#include "nusystematics/systproviders/MINERvAq0q3Weighting_tool.hh"
#include "nusystematics/systproviders/MKSinglePiTemplate_tool.hh"
#include "nusystematics/systproviders/MiscInteractionSysts_tool.hh"
#include "nusystematics/systproviders/NOvAStyleNonResPionNorm_tool.hh"
#include "nusystematics/systproviders/CCQERPAReweight_tool.hh"
#include "nusystematics/systproviders/DIRTExample_tool.hh"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace nusyst {

NEW_SYSTTOOLS_EXCEPT(unknown_nusyst_systprovider);

inline std::unique_ptr<IGENIESystProvider_tool>
make_instance(fhicl::ParameterSet const &paramset) {
  std::string tool_type = paramset.get<std::string>("tool_type");

  if (tool_type == "GENIEReWeight") {
    return std::make_unique<GENIEReWeight>(paramset);
  } else if (tool_type == "MKSinglePiTemplate") {
    return std::make_unique<MKSinglePiTemplate>(paramset);
  } else if (tool_type == "MINERvAq0q3Weighting") {
    return std::make_unique<MINERvAq0q3Weighting>(paramset);
  } else if (tool_type == "NOvAStyleNonResPionNorm") {
    return std::make_unique<NOvAStyleNonResPionNorm>(paramset);
  } else if (tool_type == "MiscInteractionSysts") {
    return std::make_unique<MiscInteractionSysts>(paramset);
  } else if (tool_type == "BeRPAWeight") {
    return std::make_unique<BeRPAWeight>(paramset);
  } else if (tool_type == "MINERvAE2p2h") {
    return std::make_unique<MINERvAE2p2h>(paramset);
  } else if (tool_type == "EbLepMomShift") {
    return std::make_unique<EbLepMomShift>(paramset);
  } else if (tool_type == "FSILikeEAvailSmearing") {
    return std::make_unique<FSILikeEAvailSmearing>(paramset);
  } else if (tool_type == "CCQERPAReweight") {
    return std::make_unique<CCQERPAReweight>(paramset);
  } else if (tool_type == "DIRTExample") {
    return std::make_unique<DIRTExample>(paramset);
  } else {
    throw unknown_nusyst_systprovider()
        << "[ERROR]: Unknown tool type: " << std::quoted(tool_type);
  }
}

} // namespace nusyst

#endif
