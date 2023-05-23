#ifndef nusystematics_RESPONSE_CALCULATORS_CCQERPAReweightCalculator_HH_SEEN
#define nusystematics_RESPONSE_CALCULATORS_CCQERPAReweightCalculator_HH_SEEN

#include "systematicstools/interface/types.hh"

#include "systematicstools/interpreters/PolyResponse.hh"

#include "systematicstools/utility/ROOTUtility.hh"
#include "systematicstools/utility/exceptions.hh"

#include "fhiclcpp/ParameterSet.h"

#ifndef NO_ART
#include "cetlib/search_path.h"
#endif

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSpline.h"

NEW_SYSTTOOLS_EXCEPT(invalid_CCQE_RPA_tweak);

namespace nusyst {

  class CCQERPAReweightCalculator{

    enum ENuRange {
      LowE = 0,
      HighE = 1,
    };

  protected:

    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithRPAXSec;
    std::map<int, std::unique_ptr<TH3D>> map_ENuRange_to_WithoutRPAXSec;

    double x_FirstBinCenter, x_LastBinCenter;
    double y_FirstBinCenter, y_LastBinCenter;
    double z_FirstBinCenter, z_LastBinCenter;

  public:

    CCQERPAReweightCalculator(fhicl::ParameterSet const &InputManifest) {
      LoadInputHistograms(InputManifest);
    }
    ~CCQERPAReweightCalculator(){}

    void LoadInputHistograms(fhicl::ParameterSet const &ps);

    double GetRPAReweight(double Enu_GeV, double P_GeV, double CTheta, double parameter_value);
    std::string GetCalculatorName() const { return "CCQERPAReweightCalculator"; }

  };

  inline double CCQERPAReweightCalculator::GetRPAReweight(double Enu_GeV, double P_GeV, double CTheta, double parameter_value){

    int enu_range = (Enu_GeV<2.1) ? 0 : 1;

    //printf("[CCQERPAReweightCalculator::GetRPAReweight] (Enu_GeV, P_GeV, CTheta) = (%1.3f, %1.3f, %1.3f), enu_range = %d\n", Enu_GeV, P_GeV, CTheta, enu_range);

    static double Enu_GeV_epsil = 1E-6;
    double Enu_GeV_ForInterp = Enu_GeV;
    Enu_GeV_ForInterp = std::max( Enu_GeV_ForInterp, x_FirstBinCenter + Enu_GeV_epsil );
    Enu_GeV_ForInterp = std::min( Enu_GeV_ForInterp, x_LastBinCenter - Enu_GeV_epsil );

    static double P_GeV_epsil = 1E-6;
    double P_GeV_ForInterp = P_GeV;
    P_GeV_ForInterp = std::max( P_GeV_ForInterp, y_FirstBinCenter + P_GeV_epsil );
    P_GeV_ForInterp = std::min( P_GeV_ForInterp, y_LastBinCenter - P_GeV_epsil );

    static double CTheta_epsil = 1E-6;
    double CTheta_ForInterp = CTheta;
    CTheta_ForInterp = std::max( CTheta_ForInterp, z_FirstBinCenter + CTheta_epsil );
    CTheta_ForInterp = std::min( CTheta_ForInterp, z_LastBinCenter - CTheta_epsil );

    //printf("[CCQERPAReweightCalculator::GetRPAReweight] -> (Enu_GeV, P_GeV, CTheta) = (%1.3f, %1.3f, %1.3f)\n", Enu_GeV_ForInterp, P_GeV_ForInterp, CTheta_ForInterp);
    double xsec_WithRPA = map_ENuRange_to_WithRPAXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, P_GeV_ForInterp, CTheta_ForInterp); // CV
    double xsec_WithoutRPA = map_ENuRange_to_WithoutRPAXSec[enu_range]->Interpolate(Enu_GeV_ForInterp, P_GeV_ForInterp, CTheta_ForInterp);

    double weight = ( xsec_WithRPA * (1.-parameter_value) + xsec_WithoutRPA * parameter_value ) / xsec_WithRPA;

    //std::cout << "[CCQERPAReweightCalculator] weight = " << weight << std::endl;

    return weight;

  }

  inline void CCQERPAReweightCalculator::LoadInputHistograms(fhicl::ParameterSet const &ps) {

#ifndef NO_ART
    bool use_stashcache = ps.get<bool>("use_FW_SEARCH_PATH", false);
#endif

    std::string const &default_root_file = ps.get<std::string>("input_file", "");

    for (fhicl::ParameterSet const &val_config :
         ps.get<std::vector<fhicl::ParameterSet>>("inputs")) {
      std::string hName = val_config.get<std::string>("name");
      std::string input_hist = val_config.get<std::string>("input_hist");
      std::string input_file = val_config.get<std::string>("input_file", default_root_file);

#ifndef NO_ART

      if (use_stashcache) {
        std::string stashcache_file;
        cet::search_path sp("FW_SEARCH_PATH");
        if (!sp.find_file(input_file, stashcache_file)) {
          char *fw = getenv("FW_SEARCH_PATH");
          std::string fw_str("");
          if (fw) {
            fw_str = fw;
          }
          throw invalid_tfile()
              << "[ERROR]: Failed to find file: " << input_file
              << ", on stashcache. (FW_SEARCH_PATH=\"" << fw_str << "\")";
        }
        input_file = stashcache_file;
      }

#endif

      if(hName=="LowE_WithRPA"){
        map_ENuRange_to_WithRPAXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="LowE_WithoutRPA"){
        map_ENuRange_to_WithoutRPAXSec[0] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="HighE_WithRPA"){
        map_ENuRange_to_WithRPAXSec[1] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }
      else if(hName=="HighE_WithoutRPA"){
        map_ENuRange_to_WithoutRPAXSec[1] = std::unique_ptr<TH3D>( GetHistogram<TH3D>(input_file, input_hist) );
      }

    }

    const auto& h_WithRPAXsec = map_ENuRange_to_WithRPAXSec[0];

    const auto& XAxis_WithRPAXsec = h_WithRPAXsec->GetXaxis();
    const auto& YAxis_WithRPAXsec = h_WithRPAXsec->GetYaxis();
    const auto& ZAxis_WithRPAXsec = h_WithRPAXsec->GetZaxis();

    x_FirstBinCenter = XAxis_WithRPAXsec->GetBinCenter(1);
    x_LastBinCenter = XAxis_WithRPAXsec->GetBinCenter( XAxis_WithRPAXsec->GetNbins() );

    y_FirstBinCenter = YAxis_WithRPAXsec->GetBinCenter(1);
    y_LastBinCenter = YAxis_WithRPAXsec->GetBinCenter( YAxis_WithRPAXsec->GetNbins() );

    z_FirstBinCenter = ZAxis_WithRPAXsec->GetBinCenter(1);
    z_LastBinCenter = ZAxis_WithRPAXsec->GetBinCenter( ZAxis_WithRPAXsec->GetNbins() );

  }


} // namespace nusyst

#endif
