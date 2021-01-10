#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/TRSM.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <string>

namespace ScannerS {
namespace Interfaces {
namespace HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInputEffC;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace HiggsBoundsSignals
} // namespace Interfaces

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief Broken phase of the TRSM
 *
 * The implementation follows the conventions of
 * [1908.08554](https://arxiv.org/abs/1908.08554)
 */
class TRSMBroken : public TRSM {
public:
  //! short description
  static constexpr auto description = "TRSM broken phase";
  //! neutral Higgs names, \f$ h_1,h_2,h_3 \f$
  static constexpr std::array namesHzero{"H1", "H2", "H3"};

  /**
   * @brief Input parametrization in terms of mixing angles.
   *
   * No ordering of \f$M_a, M_b, M_c\f$ is assumed. Instead, the input mixing
   * matrix \f$R^\mathrm{in}\f$ in the basis \f$(h_a, h_b, h_c)\f$ is calculated
   * using the input mixing angles \f$\theta^i_1, \theta^i_2, \theta^i_3\f$. The
   * mass ordered states \f$h_{1,2,3}\f$ and the mass ordered mixing matrix
   * \f$R\f$ with the corresponding mixing angles are automaticcally obtained
   * from this input.
   *
   * The mixing matrix parametrization in
   * [1908.08554](https://arxiv.org/abs/1908.08554) differs from the one in
   * Utilities::MixMat3d by a minus sign for every mixing angle, ie the input
   * mixing matrix is generated using `Utilities::MixMat3d(-t1,-t2,-t3)`.
   */
  struct AngleInput {
    double mHa; //!< \f$ M_a \f$, neutral Higgs mass in GeV
    double mHb; //!< \f$ M_b \f$, neutral Higgs mass in GeV
    double mHc; //!< \f$ M_c \f$, neutral Higgs mass in GeV
    double t1;  //!< \f$ \theta^i_{hS} \f$, first input mixing angle
    double t2;  //!< \f$ \theta^i_{hX} \f$, second input mixing angle
    double t3;  //!< \f$ \theta^i_{SX} \f$, third input mixing angle
    double v;   //!< EW vev in GeV
    double vs;  //!< \f$ v_S \f$ in GeV
    double vx;  //!< \f$ v_X \f$ in GeV
  };

  //! A broken-phase TRSM parameter point
  struct ParameterPoint {
    //! mass-ordered neutral Higgs masses \f$ M_1, M_2, M_3 \f$ in GeV
    const std::array<double, nHzero> mHi;
    //! mass-ordered Higgs mixing matrix \f$R\f$
    const Eigen::Matrix3d R;
    //! corresponding mixing angles \f$ \theta_{hS},\theta_{hX},\theta_{SX}\f$
    const std::array<double, 3> theta;
    //! EW vev in GeV
    const double v;
    //! \f$ v_S \f$ in GeV
    const double vs;
    //! \f$ v_X \f$ in GeV
    const double vx;
    //! quartic potential parameters \f$ \lambda_\Phi, \lambda_S,
    //! \lambda_X, \lambda_{\Phi S}, \lambda_{\Phi X}, \lambda_{SX} \f$
    const std::array<double, 6> L;
    //! doublet mass term \f$\mu_\Phi^2\f$ in \f$\mathrm{GeV}^2\f$
    const double muHsq;
    //! singlet mass term \f$\mu_S^2\f$ in \f$\mathrm{GeV}^2\f$
    const double muSsq;
    //! singlet mass term \f$\mu_X^2\f$ in \f$\mathrm{GeV}^2\f$
    const double muXsq;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    explicit ParameterPoint(const AngleInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1",     "mH2", "mH3", "R11",   "R12",   "R13",     "R21",
        "R22",     "R23", "R31", "R32",   "R33",   "thetahS", "thetahX",
        "thetaSX", "v",   "vs",  "vx",    "LH",    "LS",      "LX",
        "LHS",     "LHX", "LSX", "muHsq", "muSsq", "muXsq"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Model implementation for Constraints::STU
   *
   * Uses the TRSM::STUInput implementation.
   *
   * @param p the parameter point
   * @return Constraints::STUDetail::STUParameters input parameters for the STU
   * calculation
   */
  static Constraints::STUDetail::STUParameters
  STUInput(const ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * Uses the effective coupling approximation for the HiggsBounds input. Calls
   * CalculateBRs and CalculateCXNs, which utilise the SM predictions tabulated
   * in HiggsBounds and store their results.
   *
   * @param p the parameter point
   * @param hbhs the HiggsBoundsSignals object to access tabulated values
   * @return Interfaces::HiggsBoundsSignals::HBInput<nHzero, nHplus> input for
   * HiggsBounds
   */
  static Interfaces::HiggsBoundsSignals::HBInputEffC<nHzero, nHplus>
  HiggsBoundsInput(
      ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &hbhs);

  /**
   * @brief Calculate the Higgs branching ratios.
   *
   * Uses the tabulated SM-like BRs from HiggsBounds and the tree-level decay
   * widths for Higgs-to-Higgs decays implemented in Tools::ScalarWidths. Calls
   * TripleHCoups and stores all triple-Higgs couplings as `c_HaHbHc`. Also
   * stores all branching ratios of the scalars as `b_Ha_XX`, where `XX` denotes
   * a pair of SM particles or a pair of `HbHc` (only if the decay is possible
   * by the mass-ordering) as well as the total widths as `w_Ha`. All of this
   * for `Ha`, `Hb`, and `Hp` in #namesHzero.
   *
   * @param p the parameter point
   * @param hbhs the HiggsBoundsSignals object to access tabulated values
   */
  static void CalculateBRs(
      ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &hbhs);

  /**
   * @brief Calculates 13 TeV LHC production cross sections.
   *
   * Rescales the SM-like cross sections tabulated in HiggsBounds and stores
   * them named as in Interfaces::HiggsBoundsSignals::SMCxn with `h` replaced by
   * `Ha` for `Ha` in #namesHzero.
   *
   * @param p the parameter point
   * @param hbhs the HiggsBoundsSignals object to access tabulated values
   */
  static void CalculateCXNs(
      ParameterPoint &p,
      const Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<nHzero, nHplus>
          &hbhs);
};
} // namespace Models
} // namespace ScannerS
