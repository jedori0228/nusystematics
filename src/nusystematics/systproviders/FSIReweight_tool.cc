#include "nusystematics/systproviders/FSIReweight_tool.hh"

#include "nusystematics/utility/exceptions.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "Framework/GHEP/GHepParticle.h"

#include "TLorentzVector.h"

using namespace systtools;
using namespace nusyst;
using namespace fhicl;

FSIReweight::FSIReweight(ParameterSet const &params)
    : IGENIESystProvider_tool(params),
      fsiReweightCalculator(nullptr),
      ResponseParameterIdx(systtools::kParamUnhandled<size_t>),
      valid_file(nullptr), valid_tree(nullptr) {}

SystMetaData FSIReweight::BuildSystMetaData(ParameterSet const &cfg,
                                                     paramId_t firstId) {

  std::cout << "[FSIReweight::BuildSystMetaData] called" << std::endl;

  SystMetaData smd;

  systtools::SystParamHeader phdr;
  if (ParseFhiclToolConfigurationParameter(cfg, "FSIReweight",
                                                 phdr, firstId)) {
    phdr.systParamId = firstId++;
    smd.push_back(phdr);
  }

  fhicl::ParameterSet templateManifest =
      cfg.get<fhicl::ParameterSet>("FSIReweight_input_manifest");
  tool_options.put("FSIReweight_input_manifest", templateManifest);

  // OPTION_IN_CONF_FILE can be defined in the configuration file
  // then it is copied to tool_option when running "GenerateSystProviderConfig" to generation paramHeader

  fill_valid_tree = cfg.get<bool>("fill_valid_tree", false);
  tool_options.put("fill_valid_tree", fill_valid_tree);

  return smd;
}

bool FSIReweight::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {

  std::cout << "[FSIReweight::SetupResponseCalculator] called" << std::endl;

  fhicl::ParameterSet templateManifest =
      tool_options.get<fhicl::ParameterSet>("FSIReweight_input_manifest");

  std::string kin_Y_str(""), kin_Z_str("");
  kin_Y_str = "q3";
  kin_Z_str = "q0";

  std::cout << "[FSIReweight::SetupResponseCalculator] Template binnings are" << std::endl;
  std::cout << "[FSIReweight::SetupResponseCalculator] x: E_nu" << std::endl;
  std::cout << "[FSIReweight::SetupResponseCalculator] y: " << kin_Y_str << std::endl;
  std::cout << "[FSIReweight::SetupResponseCalculator] z: " << kin_Z_str << std::endl;

  fsiReweightCalculator = std::make_unique<FSIReweightCalculator>( templateManifest );

  ResponseParameterIdx =
      GetParamIndex(GetSystMetaData(), "FSIReweight");

  fill_valid_tree = tool_options.get<bool>("fill_valid_tree", false);
  if (fill_valid_tree) {
    InitValidTree();
  }

  return true;
}

event_unit_response_t
FSIReweight::GetEventResponse(genie::EventRecord const &ev) {

  // when the event is not applicable for this type of reweighting,
  // use GetDefaultEventResponse() to return an auto-1.-filled vector

  //if (!ev.Summary()->ProcInfo().IsQuasiElastic() ||
  //    !ev.Summary()->ProcInfo().IsWeakCC()) {
  //  return this->GetDefaultEventResponse();
  //}

  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  genie::GHepParticle *ISNuc = ev.TargetNucleus();
  genie::GHepParticle *ISHitNucleon = ev.HitNucleon();

  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }

  //std::cout<<ev.Summary()->AsString()<<std::flush;
  //ev.Print();
  //std::cout<<std::endl;

  vector<genie::GHepParticle*> priorFShadron_list;
  /*for (int ip=ISHitNucleon->FirstDaughter(); ip<=ISHitNucleon->LastDaughter(); ip++) {
    genie::GHepParticle *priorFShadron = ev.Particle(ip);
    cout<<priorFShadron->Pdg()<<"\t"<<priorFShadron->Status()<<endl;
    if (priorFShadron->Status() == 14)
      priorFShadron_list.push_back(priorFShadron);
  }*/
  TObjArrayIter piter(&ev);
  genie::GHepParticle * par = 0;
  while( (par = (genie::GHepParticle *) piter.Next()) ) {
    if (par->Status() == 14) {
      priorFShadron_list.push_back(par);
    }
  } // get all prior FS hadrons
  
  TLorentzVector FSLepP4 = *FSLep->P4(); // l
  TLorentzVector ISLepP4 = *ISLep->P4(); // nu
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);

  double AngleLeps = FSLepP4.Vect().Angle( ISLepP4.Vect() );
  double CAngleLeps = TMath::Cos(AngleLeps);
  double SAngleLeps = TMath::Sin(AngleLeps);
  if(AngleLeps>=M_PI/2.) SAngleLeps *= -1.;

  std::array<double, 2> bin_kin;
  bin_kin = {emTransfer.Vect().Mag(), emTransfer.E()};

  // now make the output
  systtools::event_unit_response_t resp;

  SystParamHeader const &hdr = GetSystMetaData()[ResponseParameterIdx];

  resp.push_back( {hdr.systParamId, {}} );
  for (double var : hdr.paramVariations) {
    double weight = 1;
    for (auto IShad : priorFShadron_list) {
      if (IShad->Pdg() != 2212) { // only for proton for test
        //weight = 1;
        //break;
        continue;
      }
      vector<genie::GHepParticle*> FSdaughter_list;
      IShad->Status();
      get_FS_daughters(IShad, FSdaughter_list, ev);
      double KEini = IShad->KinE();
      double Ehad = 0;
      for (auto FSpar : FSdaughter_list) {
        int pdg = FSpar->Pdg();
        double totE = FSpar->E();
        double kinE = FSpar->KinE();
        if ( pdg==211 || pdg==-211 || pdg==111) // pion
          Ehad += totE;
        else if ( pdg==2212 ) // proton
          Ehad += kinE;
        else if ( pdg==2112 ) // neutron
          ;
        else if ( pdg==22 ) // gamma
          Ehad += totE;
        else cout<<"$#@ "<<FSpar->Name()<<endl;
      }
      double Ebias = (KEini - Ehad) / KEini;
      //cout<<IShad->Name()<<": KEini "<<KEini<<"; Ehad "<<Ehad<<"; Ebias "<<Ebias<<endl;
      double this_reweight = 1;//fsiReweightCalculator->GetFSIReweight(KEini, Ebias, var, IShad->Pdg());
      weight *= this_reweight;
    }
    resp.back().responses.push_back( weight );
  }

  if (fill_valid_tree) {

    pdgfslep = ev.FinalStatePrimaryLepton()->Pdg();
    momfslep = FSLepP4.Vect().Mag();
    cthetafslep = FSLepP4.Vect().CosTheta();

    Pdgnu = ISLep->Pdg();
    NEUTMode = 0;
    if (ev.Summary()->ProcInfo().IsMEC() &&
        ev.Summary()->ProcInfo().IsWeakCC()) {
      NEUTMode = (Pdgnu > 0) ? 2 : -2;
    } else {
      NEUTMode = genie::utils::ghep::NeutReactionCode(&ev);
    }

    QELikeTarget_t qel_targ = GetQELikeTarget(ev);
    QELTarget = e2i(qel_targ);

    Enu = ISLepP4.E();
    Q2 = -emTransfer.Mag2();
    W = ev.Summary()->Kine().W(true);
    q0 = emTransfer.E();
    q3 = emTransfer.Vect().Mag();

    valid_tree->Fill();
  }

  return resp;

}

std::string FSIReweight::AsString() { return ""; }

void FSIReweight::InitValidTree() {

  valid_file = new TFile("MINERvAq3q0WeightValid.root", "RECREATE");
  valid_tree = new TTree("valid_tree", "");

  valid_tree->Branch("NEUTMode", &NEUTMode);
  valid_tree->Branch("QELTarget", &QELTarget);
  valid_tree->Branch("Enu", &Enu);
  valid_tree->Branch("Pdg_nu", &Pdgnu);
  valid_tree->Branch("Pdg_FSLep", &pdgfslep);
  valid_tree->Branch("P_FSLep", &momfslep);
  valid_tree->Branch("CosTheta_FSLep", &cthetafslep);
  valid_tree->Branch("Q2", &Q2);
  valid_tree->Branch("W", &W);
  valid_tree->Branch("q0", &q0);
  valid_tree->Branch("q3", &q3);
}

FSIReweight::~FSIReweight() {
  if (valid_file) {
    valid_tree->SetDirectory(valid_file);
    valid_file->Write();
    valid_file->Close();
    delete valid_file;
  }
}
