#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"
#include "ScannerS/Models/N2HDM.hpp"
#include "ScannerS/Models/TwoHDM.hpp"
#include <Eigen/Core>
#include <array>
#include <map>
#include <string>
#include <vector>

namespace ScannerS {
namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Tools {
class SushiTables;
} // namespace Tools

namespace Models {
/**
 * @brief The dark doublet phase of the N2HDM
 *
 * Conventions based on [1805.00966](https://arxiv.org/abs/1805.00966).
 */
class N2HDMDarkS : public N2HDM {
public:
  //! short description
  static constexpr auto description = "N2HDM dark singlet phase";
  //! names of the neutral Higgs bosons, \f$ A, H_1, H_2, H_D \f$
  static constexpr std::array namesHzero{"A", "H1", "H2", "HD"};
  //! charged Higgs names, \f$ H^\pm \f$
  static constexpr std::array namesHplus{"Hp"};
  //! number of visible-sector neutral Higgs bosons
  static constexpr int nHzeroVisible = nHzero - 1;
  //! Yukawa types
  using Yuk = TwoHDM::Yuk;

  /**
   * @brief Input parametrization in terms of mixing angles
   *
   * No ordering of \f$m_{H_a},m_{H_b}\f$ is assumed. Instead, the input mixing
   * matrix \f$R^\mathrm{in}\f$ in the basis \f$(H_a, H_b, H_D)\f$ is calculated using the
   * input mixing angle\f$\alpha^\mathrm{in}\f$. The mass ordered states \f$H_1\f$ and
   * \f$H_2\f$ and the corresponding mixing matrix \f$R\f$ and mixing angle
   * \f$alpha\f$ are automatically obtained.
   */
  struct AngleInput {
    double mHa;   //!< \f$ m_{H_a} \f$, neutral CP-even Higgs mass in GeV
    double mHb;   //!< \f$ m_{H_b} \f$, neutral CP-even Higgs mass in GeV
    double mA;    //!< \f$ m_A \f$, neutral CP-odd Higgs mass in GeV
    double mHp;   //!< \f$ m_{H^\pm} \f$, charged Higgs mass in GeV
    double mHD;   //!< \f$ m_{H_D} \f$, dark singlet Higgs mass in GeV
    double tbeta; //!< \f$\tan\beta\f$, ratio of vacuum expectation values
    double alpha; //!< input mixing angle \f$ \alpha^\mathrm{in} \f$
    double m12sq; //!< \f$ m_{12}^2\f$ in \f$\mathrm{GeV}^2\f$
    double L6;    //!< \f$\lambda_6\f$
    double L7;    //!< \f$\lambda_7\f$
    double L8;    //!< \f$\lambda_8\f$
    Yuk type;     //!< the Yukawa type
    double v;     //!< EW vev in GeV
  };

  //! Parameter point of the dark singlet phase N2HDM
  struct ParameterPoint {
    //! ordered visible neutral CP-even Higgs masses \f$m_{H_{1,2}}\f$ in GeV
    const std::array<double, 2> mHi;
    //! CP-odd Higgs mass \f$m_A\f$ in GeV
    const double mA;
    //! charged Higgs mass \f$m_{H^\pm}\f$ in GeV
    const double mHp;
    //! dark singlet Higgs mass \f$ m_{H_D} \f$ in GeV
    const double mHD;
    //! \f$\tan\beta\f$, ratio of doublet vacuum expectation values
    const double tbeta;
    //! CP-even neutral visible sector mixing angle \f$\alpha\f$
    const double alpha;
    //! mass-ordered (for \f$H_{1,2}\f$) mixing matrix \f$R\f$ (with \f$H_D\f$)
    const Eigen::Matrix3d R;
    //! Yukawa type
    const Yuk type;
    //! EW vev in GeV
    const double v;
    //! soft \f$Z_2\f$ breaking parameter \f$m_{12}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m12sq;
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
        "mH1",     "mH2", "mA",    "mHp",   "mHD",   "tbeta", "alpha", "R11",
        "R12",     "R13", "R21",   "R22",   "R23",   "R31",   "R32",   "R33",
        "yuktype", "v",   "m12sq", "L1",    "L2",    "L3",    "L4",    "L5",
        "L6",      "L7",  "L8",    "m11sq", "m22sq", "mssq"};

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
   * p.data with the names given in AnyHdecay::Hdecay::n2hdmDarkSingletKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculates Higgs couplings.
   *
   * Calculates and stores the neutral Higgs-gauge (\f$c(H_iVV)\f$ as `c_HiVV`
   * and \f$c(H_iAZ)\f$ as `c_HiAZ`) and Higgs-fermion (\f$c(H_if\bar{f})\f$ as
   * `c_Hiuu_e`, `c_Hidd_e`, `c_Hill_e` and the CP-odd \f$c(Af\bar{f})\f$ as
   * `c_Auu_o`, `c_Add_o`, `c_All_o`) couplings for `Hi` in `{"H1","H2"}`.
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
   * \f$pp\to b\bar{b}H_i\f$ cross sections as `x_Hi_bbH` for `Hi` in
   * #namesHzero (except for `HD`). Uses Tools::SushiTables.
   *
   * @param p the parameter point
   */
  static void CalcCXNs(ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * **Requires RunHdecay and CalcCouplings to be called beforehand.**
   *
   * Uses the HiggsBounds tables to obtain and store the 13TeV LHC \f$pp\to
   * tbH^\pm\f$ cross section as `x_tHpm` (in pb).
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
   *
   * **The MicrOMEGAs implementation only support Yukawa type I**
   *
   * @param p the parameter point
   * @return std::map<std::string, double> of MicrOMEGAs input parameters
   */
  static std::map<std::string, double> MOInput(const ParameterPoint &p);
  //! Constraints::DM MicrOMEGAs model name
  static constexpr auto micromegasModelName = "N2HDMDarkS_T1";

private:
  static const Tools::SushiTables cxnH0_;
};
} // namespace Models
} // namespace ScannerS
