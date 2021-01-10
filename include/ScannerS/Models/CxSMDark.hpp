#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/CxSM.hpp"
#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace ScannerS {
namespace Interfaces {
namespace HiggsBoundsSignals {
template <int nHneut, int nHplus> struct HBInput;
template <size_t nHneut, size_t nHplus> class HiggsBoundsSignals;
} // namespace HiggsBoundsSignals
} // namespace Interfaces

namespace Tools {
class SushiTables;
}

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief The dark phase of the CxSM
 *
 * Implemented by Tizian Römer. The implementation follows the conventions of
 * [1512.05355](https://arxiv.org/abs/1512.05355) but assumes \f$a_1=0\f$ and
 * the dark Higgs is called \f$h_X\f$ instead of \f$h_A\f$ to avoid confusing it
 * with a pseudoscalar.
 */
class CxSMDark : public CxSM {
public:
  //! short description
  static constexpr auto description = "CxSM dark phase";
  //! neutral Higgs names, \f$ h_1,h_2,h_X\f$
  static constexpr std::array namesHzero{"H1", "H2", "HX"};
  //! number of visible-sector neutral Higgs bosons
  static constexpr int nHzeroVisible = nHzero - 1;

  /**
   * @brief Input parametrization in terms of a mixing angle.
   *
   * No ordering of \f$m_{a,b}\f$ is required. The mass ordered mixing angle
   * \f$\alpha\f$ is automatically obtained from the input mixing angle
   * \f$\alpha^\mathrm{in}\f$ in the basis \f$(h_a, h_b)\f$.
   */
  struct AngleInput {
    double mHa;   //!< \f$ m_a \f$, neutral Higgs mass in GeV
    double mHb;   //!< \f$ m_b \f$, neutral Higgs mass in GeV
    double mHX;   //!< \f$ m_X \f$, dark Higgs mass in GeV
    double alpha; //!< input mixing angle \f$\alpha^\mathrm{in}\f$
    double v;     //!< EW vev in GeV
    double vs;    //!< \f$ v_S \f$ in GeV
  };

  //! A dark-phase CxSM parameter point
  struct ParameterPoint {
    //! light neutral CP-even Higgs mass \f$m_1\f$ in GeV
    const double mHl;
    //! heavy neutral CP-even Higgs mass \f$m_2\f$ in GeV
    const double mHh;
    //! dark Higgs mass \f$m_X\f$ in GeV
    const double mHX;
    //! visible sector mixing angle
    const double alpha;
    //! EW vev in GeV
    const double v;
    //! real singlet vev \f$ v_S \f$ in GeV
    const double vs;
    //! quartic potential parameters \f$ \lambda, d_2, \delta_2 \f$
    const std::array<double, 3> L;
    //! doublet mass term \f$ m^2 \f$
    const double msq;
    //! soft \f$U(1)\f$-breaking singlet mass term \f$b_1\f$
    const double b1;
    //! \f$U(1)\f$-invariant singlet mass term \f$b_2\f$
    const double b2;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    ParameterPoint(const AngleInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1",    "mH2", "mHX", "alpha", "lambda", "d2",
        "delta2", "msq", "b2",  "b1",    "v",      "vs"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Model implementation for Constraints::STU
   * @param p the parameter point
   * @return Constraints::STUDetail::STUParameters input parameters for the STU
   * calculation
   */
  static Constraints::STUDetail::STUParameters
  STUInput(const ParameterPoint &p);

  /**
   * @brief Runs Hdecay for the given parameter point.
   *
   * Uses the AnyHdecay interface to call sHDECAY. Results are stored in p.data
   * with the names given in AnyHdecay::Hdecay::cxsmDarkKeys.
   *
   * @param p The parameter point.
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculates Higgs couplings.
   *
   * Calculates the mixing matrix elements/kappa factors \f$R_{11}\f$ (`RH11`)
   * and \f$R_{21}\f$ (`RH21`).
   *
   * @param p The parameter point.
   */
  static void CalcCouplings(ParameterPoint &p);

  /**
   * @brief Calculate and store some LHC production cross sections
   *
   * Stores the 13TeV LHC \f$gg\to H_i\f$ cross sections as `x_Hi_ggH` and the
   * \f$pp\to b\bar{b}H_i\f$ cross sections as `x_Hi_bbH` for `Hi` in
   * `{"H1","H2"}`. Uses Tools::SushiTables.
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
      const ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &hbhs);

  /**
   * @brief Model implementation for Constraints::DM
   * @param p the parameter point
   * @return std::map<std::string, double> of MicrOMEGAs input
   * parameters
   */
  static std::map<std::string, double> MOInput(const ParameterPoint &p);
  //! Constraints::DM MicrOMEGAs model name
  static constexpr auto micromegasModelName = "DCxSM";

  //! BSMPT model name for Constraints::EWPT
  static constexpr auto bsmptModelName = "cxsm";

  /**
   * @brief Model implementation for Constraints::EWPT
   * @param p the parameter point
   * @return std::vector<double> input for the EWPT calculation in BSMPT
   */
  static std::vector<double> BsmptInput(const ParameterPoint &p);

private:
  static const std::array<std::string, 12> parNames_;

  static const Tools::SushiTables cxnH0_;
};
} // namespace Models
} // namespace ScannerS
