#ifndef nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_FSIReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_FSI_tweak);
NEW_SYSTTOOLS_EXCEPT(invalid_FSI_FILEPATH);
using namespace std;
namespace nusyst {

  class FSIReweightCalculator{

    enum ENuRange {
      LowE = 0,
      HighE = 1,
    };

  protected:

    TH2D *hist_hA2018_proton;
    TH2D *hist_hN2018_proton;
    TH2D *hist_INCL_proton;
    TH2D *hist_G4BC_proton;


  public:

    FSIReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~FSIReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetFSIReweight(double KEini, double Ebias, double parameter_value, int parpdg);

    std::string GetCalculatorName() const { return "FSIReweightCalculator"; }

  };

  inline double FSIReweightCalculator::GetFSIReweight(double KEini, double Ebias, double parameter_value, int parpdg){
    TH2D *hist_hA2018, *hist_hN2018, *hist_INCL, *hist_G4BC;
    if (parpdg == 2212) {
      hist_hA2018 = hist_hA2018_proton;
      hist_hN2018 = hist_hN2018_proton;
      hist_INCL = hist_INCL_proton;
      hist_G4BC = hist_G4BC_proton;
    }
    //printf("[FSIReweightCalculator::GetFSIReweight] -> (Enu_GeV, kin_Y, kin_Z) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, kin_Y_ForInterp, kin_Z_ForInterp);
    int idx_KEini = hist_hA2018->GetXaxis()->FindBin(KEini);
    int idx_Ebias = hist_hA2018->GetYaxis()->FindBin(Ebias);
    double weight_hA2018 = hist_hA2018->GetBinContent(idx_KEini, idx_Ebias); // CV
    double weight_hN2018 = hist_hN2018->GetBinContent(idx_KEini, idx_Ebias);
    //cout<<"idx_KEini "<<idx_KEini<<"; idx_Ebias "<<idx_Ebias<<endl;
    //cout<<"weight_hA2018 "<<weight_hA2018<<endl;
    //cout<<"weight_hN2018 "<<weight_hN2018<<endl;

    //printf("[FSIReweightCalculator::GetFSIReweight] xsec (With RPA, Without RPA) = (%1.3f, %1.3e)\n", weight_hA2018, weight_hN2018);

    if(weight_hA2018==0.){
      //cout<<"weight_hA2018==0."<<endl;
      return 1.;
    }

    double weight = ( weight_hA2018 * (1.-parameter_value) + weight_hN2018 * parameter_value ) / weight_hA2018;
    //cout<<"weight "<<weight<<endl;

    return weight;

  }

  inline void FSIReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

    std::string const &default_root_file = ps.get<std::string>("input_file", "");

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file); // If specified per hist, replace it

      // if it does not start with "/", find it under ${NUSYSTEMATICS_FQ_DIR}/data/
      if(input_file.find("/")!=0){
        std::string tmp_NUSYSTEMATICS_ROOT = std::getenv("nusystematics_ROOT");
        if(tmp_NUSYSTEMATICS_ROOT==""){
          throw invalid_FSI_FILEPATH() << "[ERROR]: ${nusystematics_ROOT} not set but put relative path:" << input_file;
        }
        input_file = tmp_NUSYSTEMATICS_ROOT+"/data/"+input_file;
      }

      if(hName=="hist_hA2018_proton"){
        hist_hA2018_proton = GetHistogram<TH2D>(input_file, input_hist);
      }
      else if(hName=="hist_hN2018_proton"){
        hist_hN2018_proton = GetHistogram<TH2D>(input_file, input_hist);
      }
    }
  }
} // namespace nusyst

#endif
