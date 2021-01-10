#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/CxSM.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <string>

namespace ScannerS {

namespace Interfaces::HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInput;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace Interfaces::HiggsBoundsSignals

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Tools {
class SushiTables;
}

namespace Models {
/**
 * @brief The broken phase of the CxSM
 *
 * Implemented by Philipp Basler. The implementation follows the conventions of
 * [1512.05355](https://arxiv.org/abs/1512.05355)
 */
class CxSMBroken : public CxSM {
public:
  //! short description
  static constexpr auto description = "CxSM broken phase";
  //! neutral Higgs names, \f$ h_1,h_2,h_3 \f$
  static constexpr std::array namesHzero{"H1", "H2", "H3"};

  /**
   * @brief Input parametrization in terms of mixing angles.
   *
   * The third neutral input Higgs mass \f$m_c\f$ is calculated from \f$m_a,
   * m_b\f$ and the initial mixing matrix \f$ R^\mathrm{in} \f$ in the basis
   * \f$h_a, h_b, h_c\f$ (= Utilities::MixMat3d() with the input mixing angles
   * \f$\alpha^\mathrm{in}_{1,2,3}\f$). No ordering of \f$m_{a,b,c}\f$ is
   * required. The mass-ordered mixing matrix \f$ R \f$ and the corresponding
   * mass-ordered eigenstates \f$ h_{1,2,3} \f$ and mixing angles
   * \f$\alpha_{1,2,3}\f$ are automatically obtained.
   */
  struct AngleInput {
    double mHa; //!< \f$ m_a \f$, neutral Higgs mass in GeV
    double mHb; //!< \f$ m_b \f$, neutral Higgs mass in GeV
    double a1;  //!< first mixing angle \f$\alpha^\mathrm{in}_1\f$
    double a2;  //!< second mixing angle \f$\alpha^\mathrm{in}_2\f$
    double a3;  //!< third mixing angle \f$\alpha^\mathrm{in}_3\f$
    double v;   //!< EW vev in GeV
    double vs;  //!< \f$ v_S \f$ in GeV
  };

  //! A broken-phase CxSM parameter point
  struct ParameterPoint {
    //! mass-ordered neutral Higgs masses \f$ m_1, m_2, m_3 \f$ in GeV
    const std::array<double, nHzero> mHi;
    //! mass-ordered Higgs mixing matrix \f$R\f$
    const Eigen::Matrix3d R;
    //! the corresponding mixing angles \f$\alpha_{1,2,3}\f$
    const std::array<double, 3> alpha;
    //! EW vev in GeV
    const double v;
    //! real singlet vev \f$ v_S \f$ in GeV
    const double vs;
    //! imaginary singlet vev \f$ v_A \f$ in GeV
    const double va;
    //! quartic potential parameters \f$ \lambda, d_2, \delta_2 \f$
    const std::array<double, 3> L;
    //! soft \f$U(1)\f$-breaking tadpole \f$ a_1 \f$
    const double a1;
    //! doublet mass term \f$ m^2 \f$
    const double msq;
    //! soft \f$U(1)\f$-breaking singlet mass term \f$b_1\f$
    const double b1;
    //! \f$U(1)\f$-invariant singlet mass term \f$b_2\f$
    const double b2;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    explicit ParameterPoint(const AngleInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1",    "mH2",    "mH3", "R1h",    "R1s", "R1a",    "R2h",
        "R2s",    "R2a",    "R3h", "R3s",    "R3a", "alpha1", "alpha2",
        "alpha3", "lambda", "d2",  "delta2", "msq", "b2",     "b1",
        "a1",     "v",      "vs",  "va"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Checks if the generation in the ParameterPoint constructor was
   * successful, **ALWAYS CALL THIS**.
   *
   * The construction of a ParameterPoint in the CxSM is not guaranteed to work
   * for all possible inputs. Since the third neutral Higgs mass is calculated
   * from the other two and the mixing matrix, it is not guaranteed to be
   * physical. If the input would result in a tachyonic third Higgs mass,
   * \f$m^2_z<0\f$, it is set to `-1`. Since `p.mHi` is mass-ordered, this
   * means that `p` is invalid if \f$m_1<0\f$.
   *
   * @param p the parameter point
   * @return is `p` valid
   */
  static inline bool Valid(const ParameterPoint &p) { return p.mHi[0] > 0; }

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
   * with the names given in AnyHdecay::Hdecay::cxsmBrokenKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculate and store some LHC production cross sections
   *
   * Stores the 13TeV LHC \f$gg\to H_i\f$ cross sections as `x_Hi_ggH` and the
   * \f$pp\to b\bar{b}H_i\f$ cross sections as `x_Hi_bbH` for `Hi` in
   * namesHzero. Uses Tools::SushiTables.
   *
   * @param p the parameter point
   */
  static void CalcCXNs(ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * **Requires RunHdecay to be called beforehand.**
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

  //! BSMPT model name for Constraints::EWPT
  static constexpr auto bsmptModelName = "cxsm";

  /**
   * @brief Model implementation for Constraints::EWPT
   * @param p the parameter point
   * @return std::vector<double> input for the EWPT calculation in BSMPT
   */
  static std::vector<double> BsmptInput(const ParameterPoint &p);

private:
  static const Tools::SushiTables cxnH0_;
};
} // namespace Models
} // namespace ScannerS
