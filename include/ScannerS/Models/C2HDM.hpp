#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/TwoHDM.hpp"
#include <Eigen/Core>
#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace ScannerS {
namespace Interfaces {
namespace HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInput;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace HiggsBoundsSignals
} // namespace Interfaces

namespace Tools {
class SushiTables;
} // namespace Tools

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief The CP-violating 2HDM
 *
 * This implementation follows the conventions used in
 * [1711.09419](http://arxiv.org/abs/arXiv:1711.09419).
 *
 */
class C2HDM : public TwoHDM {
public:
  //! short description
  static constexpr auto description = "CP-violating 2HDM";
  //! neutral Higgs names, \f$ H_1, H_2, H_3 \f$
  static constexpr std::array namesHzero{"H1", "H2", "H3"};

  /**
   * @brief Input parametrization in terms of mixing angles.
   *
   * The third neutral Higgs mass \f$m_{H_c}\f$ is calculated from
   * \f$m_{H_a}\f$, \f$m_{H_b}\f$, and the mixing matrix \f$ R^\mathrm{in} \f$
   * (= Utilities::MixMat3d() with the initial mixing angles \f$
   * \alpha^\mathrm{in}_{1,2,3}\f$). No mass ordering between
   * \f$m_{H_{a,b,c}}\f$ is required. The initial mixing matrix \f$
   * R^\mathrm{in} \f$ corresponds to the ordering \f$(H_a, H_b, H_c)\f$. The
   * mass-ordered mixing matrix \f$ R \f$ and the corresponding mass-ordered
   * eigenstates \f$ H_{1,2,3} \f$ and mixing angles \f$\alpha_{1,2,3}\f$ are
   * automatically obtained.
   *
   */
  struct AngleInput {
    double mHa;      //!< \f$ m_{H_a} \f$, neutral Higgs mass in GeV
    double mHb;      //!< \f$ m_{H_b} \f$, neutral Higgs mass in GeV
    double mHp;      //!< \f$ m_{H^\pm} \f$, charged Higgs mass in GeV
    double a1;       //!< first mixing angle \f$\alpha^\mathrm{in}_1\f$
    double a2;       //!< second mixing angle \f$\alpha^\mathrm{in}_2\f$
    double a3;       //!< third mixing angle \f$\alpha^\mathrm{in}_3\f$
    double tbeta;    //!< \f$\tan\beta\f$, ratio of vacuum expectation values
    double re_m12sq; //!< \f$ \Re(m_{12}^2) \f$, in \f$\mathrm{GeV}^2\f$
    Yuk type;        //!< the Yukawa type
    double v;        //!< EW vev in GeV
  };

  /**
   * @brief Physical input parametrization.
   *
   * A reparametrization using couplings instead of mixing angles to allow a
   * more straightforward selection of scan ranges compared to AngleInput when
   * considering \f$H_a\equiv h_{125}\f$.
   *
   * The AngleInput parameters \f$\alpha^\mathrm{in}_{1,2,3}\f$ are expressed
   * through \f$c^2(H_aVV)\f$, \f$c^2(H_at\bar{t})\f$, \f$R^\mathrm{in}_{y3}\f$
   * and \f$\mathrm{sign}(R^\mathrm{in}_{x3})\f$. Two signs are needed as input
   * to compensate for the *squared* input for the effective couplings. The
   * first one is \f$\mathrm{sign}(R^\mathrm{in}_{x3})\f$ and the second one is
   * fixed through the *physical* assumption, that \f$c(H_aVV)\times
   * c^e(H_at\bar{t})>0\f$, which is a known consequence of Higgs measurements
   * for \f$H_a = h_{125}\f$.
   */
  struct PhysicalInput {
    double mHa; //!< @copybrief AngleInput::mHa
    double mHb; //!< @copybrief AngleInput::mHb
    double mHp; //!< @copybrief AngleInput::mHp
    //! \f$c^2(H_aVV)\f$, squared effective coupling of \f$ H_a\f$ to W/Z bosons
    double c_HaVV_sq;
    //! \f$|c(H_at\bar{t})|^2 = (c^e_t)^2 + (c^o_t)^2\f$, squared effective
    //! coupling of \f$ H_a\f$ to top quarks
    double c_Hatt_sq;
    //! \f$\mathrm{sign}(R^\mathrm{in}_{a3})\f$, sign of the
    //! \f$R^\mathrm{in}_{a3}\f$ mixing matrix element
    int sign_Ra3;
    //! \f$R_{b3}\f$ mixing element between \f$H_b\f$ and the CP-odd field
    double Rb3;
    double tbeta;    //!< @copybrief AngleInput::tbeta
    double re_m12sq; //!< @copybrief AngleInput::re_m12sq
    Yuk type;        //!< @copybrief AngleInput::type
    double v;        //!< @copybrief AngleInput::v
  };

  //! A C2HDM parameter point
  struct ParameterPoint {
    //! mass-ordered neutral Higgs masses \f$m_{H_{1,2,3}}\f$ in GeV
    const std::array<double, nHzero> mHi;
    //! charged Higgs mass \f$m_{H^\pm}\f$ in GeV
    const double mHp;
    //! \f$\tan\beta\f$, ratio of vacuum expectation values
    const double tbeta;
    //! mass-ordered neutral mixing matrix \f$R\f$
    const Eigen::Matrix<double, nHzero, nHzero> R;
    //! the corresponding mixing angles \f$\alpha_{1,2,3}\f$
    const std::array<double, 3> alpha;
    //! quartic potential parameters
    //! \f$\lambda_{1,2,3,4},\Re(\lambda_5),\Im(\lambda_5)\f$
    const std::array<double, 6> L;
    //! soft \f$Z_2\f$ breaking parameter \f$m_{12}^2\f$ in \f$\mathrm{GeV}^2\f$
    const std::complex<double> m12sq;
    //! doublet mass term \f$m_{11}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m11sq;
    //! doublet mass term \f$m_{22}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m22sq;
    //! Yukawa type
    const Yuk type;
    //! EW vev in GeV
    const double v;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    explicit ParameterPoint(const AngleInput &in);
    //! Construct a ParameterPoint from PhysicalInput
    explicit ParameterPoint(const PhysicalInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1",     "mH2", "mH3", "mHp", "tbeta", "m12sqr", "m12sqi", "L1",
        "L2",      "L3",  "L4",  "L5r", "L5i",   "m11sq",  "m22sq",  "R11",
        "R12",     "R13", "R21", "R22", "R23",   "R31",    "R32",    "R33",
        "yuktype", "v",   "a1",  "a2",  "a3"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Checks if the generation in the ParameterPoint constructor was
   * successful, **ALWAYS CALL THIS**.
   *
   * The construction of a ParameterPoint in the C2HDM is not guaranteed to work
   * for all possible inputs. This can happen for two reasons:
   *  - Since the third neutral Higgs mass is calculated, from the other two and
   *    the mixing matrix, it is not guaranteed to be physical. If the input
   *    would result in a tachyonic third Higgs mass, \f$m^2_{H_c}<0\f$, it is
   *    set to `-1`. Since `p.mHi` is mass-ordered, this means that `p` is
   *    invalid if \f$m_{H_1}<0\f$.
   *  - When using the PhysicalInput, there are combinations of input couplings
   *    for which no solution for \f$\alpha^\mathrm{in}_{1,2,3}\f$ exists. In
   * this case, the mixing matrix is set to a constant value \f$> 1\f$ such that
   * `p` is invalid if \f$R_{11}>1\f$.
   *
   * @param p the parameter point
   * @return is `p` valid
   */
  static inline bool Valid(const ParameterPoint &p) {
    return (p.mHi[0] > 0) && (p.R(0, 0) < 1);
  }

  /**
   * @brief Model implementation for Constraints::BFB
   *
   * Uses the TwoHDM::BFB implementation.
   *
   * @param L the quartic parameters of the scalar potential
   * @return is the scalar potential at `p` bounded from below?
   */
  static bool BFB(const std::array<double, 6> &L);

  /**
   * @brief Model implementation for Constraints::Unitarity
   *
   * Uses the TwoHDM::MaxUnitarityEV implementation.
   *
   * @param L Quartic couplings
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(const std::array<double, 6> &L);

  /**
   * @brief Model implementation for Constraints::AbsoluteStability
   *
   * Implements eq (23) from [1507.05100](https://arxiv.org/abs/1507.05100).
   *
   * @param p the parameter point
   * @return is the EW vacuum at p absolutely stable?
   */
  static bool AbsoluteStability(const ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::STU
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
   * Uses the AnyHdecay interface to call c2hdm_hdecay. Results are stored in
   * p.data with the names given in AnyHdecay::Hdecay::c2hdmKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculates Higgs couplings.
   *
   * Calculates and stores the effective gauge couplings \f$c(H_iVV)\f$
   * (`c_HiVV`) and \f$c(H_iH_jZ)\f$ (`c_HiHjZ`) as well as the CP-even
   * (\f$c^e(H_if\bar{f})\f$, `_e`) and CP-odd (\f$c^o(H_if\bar{f})\f$, `_o`)
   * fermion couplings `c_Hiuu_e/o` `c_Hidd_e/o`, `c_Hill_e/o` for each `Hi` and
   * `Hj` in #namesHzero. Also calculates and stores the dimensionless (ie
   * normalized to the EW vev) neutral Higgs - charged Higgs couplings as
   * `c_HpHmHj`.
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
   * #namesHzero. Uses Tools::SushiTables.
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
   * \bar{t}bH^+\f$ cross section as `x_tHpm` (in pb).
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
   * @brief Model implementation for Constraints::ElectronEDM
   *
   * Calculates the electric dipole moment of the electron using Tools::C2HEDM.
   * Stores the individual contributions (all in \f$e\,\mathrm{cm}\f$):
   *  - the contribution of fermion loops as `edm_e_f`
   *  - the contribution of charged Higgs loops as `edm_e_Hp`
   *  - the contribution of W-boson loops as `edm_e_W`
   *  - the \f$H^\pm W^\mp\gamma\f$ contribution as `edm_e_HpW`
   *
   * The total electron EDM is stored by the constraint.
   *
   * @param p the parameter point
   * @return double prediction for the electron EDM in \f$e\,\mathrm{cm}\f$
   */
  static double CalcElectronEDM(ParameterPoint &p);

  //! BSMPT model name for Constraints::EWPT
  static constexpr auto bsmptModelName = "c2hdm";

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
