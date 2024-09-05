#pragma once

#include "nusystematics/utility/exceptions.hh"
#include "nusystematics/utility/simbUtility.hh"

#include "Framework/EventGen/EventRecord.h"

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/ParticleData/PDGUtils.h"

#include "Framework/Interaction/SppChannel.h"
#include "Framework/Interaction/ProcessInfo.h"

#include <sstream>

namespace nusyst {
    // *****************************
// TH: Taken out of MaCh3 Structs.h
// Get the mass of a particle from the PDG in GeV
inline double GetMassFromPDG(int PDG) {

    switch (abs(PDG)) {
      
    case 11:
      return 0.511E-3;
      break;
    case 13:
      return 105.658E-3;
      break;
    case 15:
      return 1.77682;
      break;
    case 22:
      return 0.;
      break;
    case 211:
      return 139.57E-3;
      break;
    case 111:
      return 134.98E-3;
      break;
    case 2112:
      return 939.565E-3;
      break;
    case 2212:
      return 938.27E-3;
      break;
    //Oxygen nucleus
    case 1000080160:
      return 14.89926;
      break;
	//eta
	case 221:
	  return 547.862E-3;
	  break;
	  //K^0 (s or l)
	case 311:
	case 130:
	case 310:
	  return 497.611E-3;
	  break;
	case 321:
	  return 493.677E-3;
	  break;
	// Lamda baryon
	case 3122:
	  return 1115.683E-3;
	  break;
    case 12:
    case 14:
    case 16:
      return 0.0;
      break;
    default:
      std::cerr << "Haven't got a saved mass for PDG " << PDG << std::endl;
      std::cerr << "Please implement me! " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    } // End switch
    
    std::cerr << "Warning, didn't catch a saved mass" << std::endl;
    return 0;
}

inline double Getq0(genie::EventRecord const &ev){
  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton();
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector FSLepP4 = *FSLep->P4();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TLorentzVector emTransfer = (ISLepP4 - FSLepP4);
  double q0 = emTransfer[3];

  return q0;
}

// TH: Adapted from NUISANCE
inline float GetPmiss(genie::EventRecord const &ev, bool preFSI) {
  // pmiss_vect is the vector difference between the neutrino momentum and the sum of final state particles momenta
  // initialize to neutrino momentum
  genie::GHepParticle *ISLep = ev.Probe();
  TLorentzVector ISLepP4 = *ISLep->P4();
  TVector3 pmiss_vect = ISLepP4.Vect();
  // Get sum of momenta for all final state particles
  TVector3 Sum_of_momenta(0, 0, 0);

  TIter event_iter(&ev);
  genie::GHepParticle *p = 0;
  TVector3 dummy(0,0,0);

  // TH: NUISANCE loop starts with i = 3
  for (int i = 3; i < ev.GetEntries(); i++){
    
    p = ev.Particle(i);

    if (preFSI){
      // TH: code from IsPrimary function in NUISANCE used to pick out primaries
      if (IsPrimary(ev, p) == false)
        continue;
    }
    else { // post-FSI loop
      // TH: Only select final state particles
      if (p->Status() != genie::kIStStableFinalState)
        continue;
    }
      
    // skip nuclear remnant
    if (abs(p->Pdg()) > 10000)
      continue;

    dummy.SetXYZ(p->Px(), p->Py(), p->Pz());
    Sum_of_momenta += dummy;

    dummy.SetXYZ(0,0,0);
  }

  pmiss_vect -= Sum_of_momenta;
  float pmiss = pmiss_vect.Mag();

  return pmiss;
}

// TH: Adapted from NUISANCE
inline float GetEmiss(genie::EventRecord const &ev, bool preFSI){
  double Emiss = -999;
  double pmiss;
  if (preFSI) pmiss = GetPmiss(ev, true);
  else pmiss = GetPmiss(ev, false);

  std::map<int, double> bindingEnergies;

  // TH: in GeV
  bindingEnergies.insert(std::pair(1000030060, 0.001 * 17.0)); //Li6   
  bindingEnergies.insert(std::pair(1000060120, 0.001 * 25.0)); //C12   
  bindingEnergies.insert(std::pair(1000080160, 0.001 * 27.0)); //O16   
  bindingEnergies.insert(std::pair(1000120240, 0.001 * 32.0)); //Mg24  
  bindingEnergies.insert(std::pair(1000180400, 0.001 * 29.5)); //Ar40  
  bindingEnergies.insert(std::pair(1000200400, 0.001 * 28.0)); //Ca40  
  bindingEnergies.insert(std::pair(1000220480, 0.001 * 30.0)); //Ti48    
  bindingEnergies.insert(std::pair(1000260560, 0.001 * 36.0)); //Fe56  
  bindingEnergies.insert(std::pair(1000280580, 0.001 * 36.0)); //Ni58  
  bindingEnergies.insert(std::pair(1000822080, 0.001 * 44.0)); //Pb208 

  int n_tgt_nucleons = ev.Summary()->InitState().Tgt().A();
  int tgt_pdg = ev.Summary()->InitState().Tgt().Pdg();
  int n_int_nucleons = 1;
  // For 2p2h subract 2 nucleons
  if (ev.Summary()->ProcInfo().IsMEC()){
    if (genie::pdg::IsNeutrino(ev.Summary()->InitState().ProbePdg()) || genie::pdg::IsAntiNeutrino(ev.Summary()->InitState().ProbePdg())){
      n_int_nucleons = 2;
    }
  }

  double M_tgt = -999;
  double M_rem = -999;
  double mass_nucleon = (GetMassFromPDG(2212) + GetMassFromPDG(2112)) * 0.5;
  if (tgt_pdg == 1000010010){
    M_tgt = mass_nucleon;
    M_rem = 0;
  } else{
    M_tgt = n_tgt_nucleons*mass_nucleon - bindingEnergies[tgt_pdg];
    M_rem = M_tgt - mass_nucleon*n_int_nucleons + (1 - (double)n_int_nucleons/n_tgt_nucleons)*bindingEnergies[tgt_pdg];
  }

  double Trem = sqrt(pmiss*pmiss + M_rem*M_rem) - M_rem;
  double Ehad = 0;

  genie::GHepParticle *p = 0;

  // TH: code adapted from FitUtils.cxx in NUISANCE
  for (int i = 3; i < ev.GetEntries(); i++){
    
    p = ev.Particle(i);

    if (preFSI){
      // Check for primary vertex if calculting pre-FSI
      if (IsPrimary(ev, p) == false)
        continue;
    }
    else { // post-FSI loop
      // TH: Only select final state particles
      if (p->Status() != genie::kIStStableFinalState)
        continue;
    }
      
    // skip nuclear remnant
    if (abs(p->Pdg()) > 10000)
      continue;
    // skip lepton
    if ( abs( ev.Particle(i)->Pdg() ) == abs(ev.Particle(0)->Pdg()) - 1 )
      continue;

    // add kinetic energy of proton of neutron
    if (p->Pdg() == 2112 || p->Pdg() == 2212){
      Ehad += p->KinE();
    } else { // add total energy of other particles
      Ehad += p->E();
    }
  }

  double q0_true = Getq0(ev);

  Emiss = q0_true - Ehad - Trem;

  return Emiss;

}

} // namespace nusyst
