#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"
#include "ScannerS/Models/N2HDM.hpp"
#include <array>
#include <map>
#include <string>
#include <vector>

namespace ScannerS {
namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {

/**
 * @brief The fully dark phase of the N2HDM
 *
 * Conventions based on the N2HDECAY
 * [manual](https://gitlab.com/api/v4/projects/jonaswittbrodt%2Fn2hdecay/jobs/artifacts/release/raw/doc/N2HDECAY.pdf?job=doc)
 *
 * @todo update with dedicated reference
 */
class N2HDMDarkSD : public N2HDM {
public:
  //! short description
  static constexpr auto description = "N2HDM fully dark phase";
  //! neutral Higgs names, \f$ H_\mathrm{SM}, H_D^D, H_D^S, A_D\f$
  static constexpr std::array namesHzero{"Hsm", "HDD", "HDS", "AD"};
  //! charged Higgs names, \f$ H^\pm_D\f$
  static constexpr std::array namesHplus{"HDp"};

  //! Input parametrization
  struct Input {
    double mHsm;  //!< \f$ m_{H_\mathrm{SM}} \f$, visible Higgs mass in GeV
    double mHDD;  //!< \f$ m_{H^D_D} \f$, dark doublet neutral Higgs mass in GeV
    double mAD;   //!< \f$ m_{A_D} \f$, dark doublet neutral Higgs mass in GeV
    double mHDp;  //!< \f$m_{H^\pm_D}\f$, dark doublet charged Higgs mass in GeV
    double mHDS;  //!< \f$ m_{H^S_D} \f$, dark singlet Higgs mass in GeV
    double m22sq; //!< \f$m_{22}^2\f$, in \f$\mathrm{GeV}^2\f$
    double mssq;  //!< \f$ m_S^2 \f$, in \f$\mathrm{GeV}^2\f$
    double L2;    //!< \f$\lambda_2\f$
    double L6;    //!< \f$\lambda_6\f$
    double L8;    //!< \f$\lambda_8\f$
    double v;     //!< EW vev in GeV
  };

  //! Parameter point of the fully dark phase N2HDM
  struct ParameterPoint {
    //! visible Higgs mass \f$ m_{H_\mathrm{SM}} \f$ in GeV
    const double mHsm;
    //! dark doublet neutral Higgs mass \f$ m_{H^D_D} \f$ in GeV
    const double mHDD;
    //! dark doublet neutral Higgs mass \f$ m_{A^D_D} \f$ in GeV
    const double mAD;
    //! dark doublet charged Higgs mass \f$ m_{H^\pm_D} \f$ in GeV
    const double mHDp;
    //! dark singlet Higgs mass \f$ m_{H^S_D} \f$ in GeV
    const double mHDS;
    //! quartic potential parameters \f$\lambda_{1,2,3,4,5,6,7,8}\f$
    const std::array<double, 8> L;
    //! doublet mass term \f$m_{11}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m11sq;
    //! doublet mass term \f$m_{22}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m22sq;
    //! singlet mass term \f$m_S^2\f$ in \f$\mathrm{GeV}^2\f$
    const double mssq;
    //! EW vev in GeV
    const double v;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from Input
    explicit ParameterPoint(const Input &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mHsm", "mHDD", "mAD", "mHDp", "mHDS",  "L1",    "L2",   "L3", "L4",
        "L5",   "L6",   "L7",  "L8",   "m11sq", "m22sq", "mssq", "v"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Model implementation for Constraints::STU
   *
   * Uses the N2HDM::STUInput implementation.
   *
   * @param p the parameter point
   * @return Constraints::STUDetail::STUParameters input parameters for the STU
   * calculation
   */
  static Constraints::STUDetail::STUParameters
  STUInput(const ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::STU
   * @param p the parameter point
   * @return is the oblique parameter approximation applicable for `p`
   */
  static bool EWPValid(const ParameterPoint &p);

  /**
   * @brief Runs Hdecay for the given parameter point.
   *
   * Uses the AnyHdecay interface to call N2HDECAY. Results are stored in
   * p.data with the names given in
   * AnyHdecay::Hdecay::n2hdmDarkSingletDoubletKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * **Requires RunHdecay to be called beforehand.**
   *
   * @param p the parameter point
   * @return Interfaces::HiggsBoundsSignals::HBInput<nHzero, nHplus> input for
   * HiggsBounds
   */
  static Interfaces::HiggsBoundsSignals::HBInput<nHzero, nHplus>
  HiggsBoundsInput(
      ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &);

  /**
   * @brief Model implementation for Constraints::VacStab
   * @param p the parameter point
   * @return std::vector<double> of EVADE input parameters
   */
  static std::vector<double> ParamsEVADE(const ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::DM
   * @param p the parameter point
   * @return std::map<std::string, double> of MicrOMEGAs input parameters
   */
  static std::map<std::string, double> MOInput(const ParameterPoint &p);
  //! Constraints::DM MicrOMEGAs model name
  static constexpr auto micromegasModelName = "N2HDMDarkSD";
};
} // namespace Models
} // namespace ScannerS
