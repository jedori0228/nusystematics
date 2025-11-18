#include "nusystematics/systproviders/MINERvASPPReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

MINERvASPPReweight::MINERvASPPReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      pidx_Q2(systtools::kParamUnhandled<size_t>),
      pidx_Tpi(systtools::kParamUnhandled<size_t>) {

}

SystMetaData MINERvASPPReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[MINERvASPPReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  for (std::string const &pname :
       {"MINERvASPP_Q2", "MINERvASPP_Tpi"}){
    systtools::SystParamHeader phdr;
    if (ParseFhiclToolConfigurationParameter(cfg, pname, phdr, firstId)) {
      phdr.systParamId = firstId++;
      smd.push_back(phdr);
    }
  }

  return smd;
}

bool MINERvASPPReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[MINERvASPPReweight::SetupResponseCalculator] called" << std::endl;

  systtools::SystMetaData const &md = GetSystMetaData();

  if (HasParam(md, "MINERvASPP_Q2")) {
    pidx_Q2 = GetParamIndex(md, "MINERvASPP_Q2");
  }


  if (HasParam(md, "MINERvASPP_Tpi")) {
    pidx_Tpi = GetParamIndex(md, "MINERvASPP_Tpi");
  }

  return true;
}

event_unit_response_t
MINERvASPPReweight::GetEventResponse(genie::EventRecord const &ev) {

  // when the event is not applicable for this type of reweighting,
  // use GetDefaultEventResponse() to return an auto-1.-filled vector

  bool IsCOH = ev.Summary()->ProcInfo().IsCoherentProduction();
  if(IsCOH){
    //std::cout << "[JSKIMDEBUG] COH" << std::endl;
    return this->GetDefaultEventResponse();
  }

  // loop over particles
  int ip=-1;
  genie::GHepParticle * p = 0;
  TIter event_iter(&ev);

  // find highest momentum final-state particle
  double max_mom_fs_pip = -999.;
  genie::GHepParticle *ptl_hm_fs_pip = 0;
  int genie_n_photons = 0;
  int genie_n_mesons = 0;
  int nPip=0;

  while ( (p = dynamic_cast<genie::GHepParticle *>(event_iter.Next())) ) {

    ip++;

    // Skip particles not rescattered by the actual hadron transport code
    int  pdgc       = p->Pdg();
    bool is_pion    = genie::pdg::IsPion   (pdgc);
    genie::GHepStatus_t ist  = p->Status();
    if( ist!=genie::kIStStableFinalState ) continue;

    // photon
    if( pdgc==22 ){
      TLorentzVector* ptl_P4 = p->P4();
      if(ptl_P4->E() > 0.010){
        genie_n_photons++;
      }
    }

    // count all mesons
    if (abs(pdgc) == 211 || //pi+-
             pdgc == 111 ||  // pi0
             abs(pdgc) == 321 || // K-
             abs(pdgc) == 323 || // K*+-
             pdgc == 130 || // KL0
             pdgc == 310 || // KS0
             pdgc == 311 || // K0
             pdgc == 313 || // K*0
             abs(pdgc) == 221 || // eta
             abs(pdgc) == 331 // eta' (958)
             ) {
      genie_n_mesons++;
    }

    // count pip and save thet momentum
    if( pdgc==genie::kPdgPiP ){

      nPip++;

      TLorentzVector* ptl_P4 = p->P4();
      double mom_mag = ptl_P4->Vect().Mag();

      if( mom_mag > max_mom_fs_pip ){
        ptl_hm_fs_pip = p;
      }

    }
    else{
      continue;
    }

  }//p

  // SPP
  if( nPip != 1 || genie_n_mesons!= 1 ){
    //std::cout << "[JSKIMDEBUG] NOT SPP" << std::endl;
    return this->GetDefaultEventResponse();
  }

  // has photon with E>10 MeV; following MINERvA CC1pip signal definition
  if(genie_n_photons!=0){
    //std::cout << "[JSKIMDEBUG] Has Photon" << std::endl;
    return this->GetDefaultEventResponse();
  }

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double this_Tpi_GeV = ptl_hm_fs_pip->KinE();
  double this_Q2_GeV2 = -emTransfer.Mag2();
  //std::cout << "Leading pi+ T = " << this_Tpi_GeV << " GeV" << std::endl;

  // now make the output
  systtools::event_unit_response_t resp;

  systtools::SystMetaData const &md = GetSystMetaData();

  int TargetA = ev.Summary()->InitState().Tgt().A();
  bool IsH = TargetA==1;

/*
    for (double var : md[pidx_SPPTpiCorrectionRW].paramVariations) {
      double this_reweight = GetSPPTpiCorrectionRW(this_Q2_GeV2, this_Tpi_GeV, var);
      if(IsH) resp.back().responses.push_back( 1. );
      else resp.back().responses.push_back( this_reweight );
    }
*/

  //std::cout << "[JSKIMDEBUG] SPP event found.." << std::endl;

  if (pidx_Q2 != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_Q2].systParamId, {}} );
    double this_Q2_RW = nusyst::MINERvASPP::GetQ2TemplateReweight(this_Q2_GeV2);
    double this_Q2_OneSigSize = this_Q2_RW-1.0;
    for (double var : md[pidx_Q2].paramVariations) {
      if(IsH) resp.back().responses.push_back( 1. );
      else{
        double this_rw = 1.0 + var * this_Q2_OneSigSize;
        resp.back().responses.push_back( this_rw );
      }
    }
  }

  if (pidx_Tpi != systtools::kParamUnhandled<size_t>) {
    resp.push_back( {md[pidx_Tpi].systParamId, {}} );
    double this_Tpi_RW = nusyst::MINERvASPP::GetTpiReweight(this_Tpi_GeV);
    double this_Tpi_OneSigSize = this_Tpi_RW - 1.0;
    for (double var : md[pidx_Tpi].paramVariations) {
      if(IsH) resp.back().responses.push_back( 1. );
      else{
        double this_rw = 1.0 + var * this_Tpi_OneSigSize;
        resp.back().responses.push_back( this_rw );
      }
    }
  }

  return resp;

}

std::string MINERvASPPReweight::AsString() { return "MINERvASPPReweight"; }

MINERvASPPReweight::~MINERvASPPReweight() {
}
