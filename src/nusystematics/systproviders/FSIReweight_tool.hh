#ifndef nusystematics_SYSTPROVIDERS_FSIReweight_TOOL_SEEN
#define nusystematics_SYSTPROVIDERS_FSIReweight_TOOL_SEEN

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "nusystematics/responsecalculators/FSIReweightCalculator.hh"
#include "nusystematics/utility/enumclass2int.hh"

#include "nusystematics/utility/GENIEUtils.hh"

#include "TFile.h"
#include "TTree.h"

#include <memory>
#include <string>

class FSIReweight : public nusyst::IGENIESystProvider_tool {

  std::unique_ptr<nusyst::FSIReweightCalculator> fsiReweightCalculator;
  TFile* savemap = new TFile("FSI_2Dmap.root", "update");
  TH2D* h_KEini_Ebias = new TH2D("h_KEini_Ebias", "Ebias vs KEini; KEini [GeV]; Ebias [GeV]", 40, 0, 2, 40, -0.9, 1.1);
public:
  explicit FSIReweight(fhicl::ParameterSet const &);

  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  std::string AsString();

  ~FSIReweight();

private:

  fhicl::ParameterSet tool_options;

  size_t ResponseParameterIdx;

  void InitValidTree();

  bool fill_valid_tree;
  TFile *valid_file;
  TTree *valid_tree;

  int NEUTMode, Pdgnu, pdgfslep, QELTarget;
  double Enu, momfslep, cthetafslep, Q2, q0, q3, W;
};

void get_FS_daughters(genie::GHepParticle* par, vector<genie::GHepParticle*>& FSdaughters, genie::EventRecord const &ev) {
  if (par->Status() == 1) {
    FSdaughters.push_back(par);
    return;
  }
  if (par->FirstDaughter()==-1) {
    cout<<"No daughter for this particle..."<<endl;
    return;
  }
  for (int ip=par->FirstDaughter(); ip<=par->LastDaughter(); ip++) {
    get_FS_daughters(ev.Particle(ip), FSdaughters, ev);
  }
}



#endif
