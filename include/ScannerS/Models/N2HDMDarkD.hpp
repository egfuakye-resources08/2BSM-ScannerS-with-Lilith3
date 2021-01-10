#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"
#include "ScannerS/Models/N2HDM.hpp"
#include <Eigen/Core>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace ScannerS {

namespace Tools {
class SushiTables;
}

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief The dark doublet phase of the N2HDM
 *
 * Conventions based on [1805.00966](https://arxiv.org/abs/1805.00966), there
 * called intert doublet phase (IDP).
 */
class N2HDMDarkD : public N2HDM {
public:
  //! short description
  static constexpr auto description = "N2HDM dark doublet phase";
  //! neutral Higgs names, \f$ H_1, H_2, H_D, A_D \f$
  static constexpr std::array namesHzero{"H1", "H2", "HD", "AD"};
  //! charged Higgs names, \f$ H^\pm_D\f$
  static constexpr std::array namesHplus{"HDp"};
  //! number of visible-sector neutral Higgs bosons
  static constexpr int nHzeroVisible = nHzero - 2;

  /**
   * @brief Input parametrization in terms of mixing angles
   *
   * No ordering of \f$m_{H_a},m_{H_b}\f$ is assumed. Instead, the input mixing
   * matrix \f$R^\mathrm{in}\f$ in the basis \f$(H_a, H_b, H_D)\f$ is calculated using the
   * input mixing angle\f$\alpha^\mathrm{in}\f$. The mass ordered states \f$H_1\f$ and
   * \f$H_2\f$ and the corresponding mixing matrix \f$ R\f$ and mixing angle
   * \f$alpha\f$ are automatically obtained.
   */
  struct AngleInput {
    double mHa;   //!< \f$ m_{H_a} \f$, neutral CP-even Higgs mass in GeV
    double mHb;   //!< \f$ m_{H_b} \f$, neutral CP-even Higgs mass in GeV
    double mHD;   //!< \f$ m_{H_D} \f$, dark doublet neutral Higgs mass in GeV
    double mAD;   //!< \f$ m_{A_D} \f$, dark doublet neutral Higgs mass in GeV
    double mHDp;  //!< \f$m_{H^\pm_D}\f$, dark doublet charged Higgs mass in GeV
    double alpha; //!< input mixing angle \f$ \alpha^\mathrm{in} \f$
    double m22sq; //!< \f$m_{22}^2\f$, in \f$\mathrm{GeV}^2\f$
    double L2;    //!< \f$\lambda_2\f$
    double L8;    //!< \f$\lambda_8\f$
    double vs;    //!< singlet vev \f$ v_s \f$ in GeV
    double v;     //!< EW vev in GeV
  };

  //! Parameter point of the dark doublet phase N2HDM
  struct ParameterPoint {
    //! mass-ordered neutral visible Higgs masses \f$m_{H_{1,2}}\f$ in GeV
    const std::array<double, nHzeroVisible> mHi;
    //! dark doublet neutral Higgs mass \f$ m_{H_D} \f$ in GeV
    const double mHD;
    //! dark doublet neutral Higgs mass \f$ m_{A_D} \f$ in GeV
    const double mAD;
    //! dark doublet charged Higgs mass \f$ m_{H^\pm_D} \f$ in GeV
    const double mHDp;
    //! visible sector mixing angle \f$\alpha\f$
    const double alpha;
    //! mass-ordered (for \f$H_{1,2}\f$) mixing matrix \f$R\f$ (with \f$H_D\f$)
    const Eigen::Matrix3d R;
    //! singlet vev \f$ v_s \f$ in GeV
    const double vs;
    //! EW vev in GeV
    const double v;
    //! quartic potential parameters \f$\lambda_{1,2,3,4,5,6,7,8}\f$
    const std::array<double, 8> L;
    //! doublet mass term \f$m_{11}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m11sq;
    //! doublet mass term \f$m_{22}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m22sq;
    //! singlet mass term \f$m_S^2\f$ in \f$\mathrm{GeV}^2\f$
    const double mssq;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    explicit ParameterPoint(const AngleInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1", "mH2", "mHD", "mAD", "mHDp",  "alpha", "R11",
        "R12", "R13", "R21", "R22", "R23",   "R31",   "R32",
        "R33", "vs",  "v",   "L1",  "L2",    "L3",    "L4",
        "L5",  "L6",  "L7",  "L8",  "m11sq", "m22sq", "mssq"};

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
   * p.data with the names given in AnyHdecay::Hdecay::n2hdmDarkDoubletKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculates Higgs couplings.
   *
   * Calculates and stores the neutral Higgs-gauge \f$c(H_iVV)\f$ (`c_HiVV`) and
   * Higgs-fermion \f$c(H_if\bar{f})\f$ (`c_Hiuu_e`, `c_Hidd_e`, `c_Hill_e`)
   * couplings for the visible Higgs bosons `Hi` in `{"H1","H2"}`.
   *
   * @param p the parameter point
   */
  static void CalcCouplings(ParameterPoint &p);

  /**
   * @brief Calculate and store some LHC production cross sections
   *
   * **Requires CalcCouplings to be called beforehand.**
   *
   * Stores the 13TeV LHC \f$gg\to H_i\f$ cross sections as `x_Hi_ggH` and the
   * \f$pp\to b\bar{b}H_i\f$ cross sections as `x_Hi_bbH` for the visible Higgs
   * bosons `Hi` in `{"H1","H2"}`. Uses Tools::SushiTables.
   *
   * @param p the parameter point
   */
  static void CalcCXNs(ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * **Requires RunHdecay and CalcCouplings to be called beforehand.**
   *
   * @param p the parameter point
   * @param hbhs the HiggsBoundsSignals object to access tabulated values
   * @return Interfaces::HiggsBoundsSignals::HBInput<nHzero, nHplus> input for
   * HiggsBounds
   */
  static Interfaces::HiggsBoundsSignals::HBInput<nHzero, nHplus>
  HiggsBoundsInput(
      ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &hbhs);

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
  static constexpr auto micromegasModelName = "N2HDMDarkD";

private:
  static const Tools::SushiTables cxnH0_;
};
} // namespace Models
} // namespace ScannerS
