#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/TwoHDM.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ScannerS {
namespace Interfaces::HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInput;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace Interfaces::HiggsBoundsSignals

namespace Tools {
class SushiTables;
} // namespace Tools

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief The real (CP-conserving) two Higgs doublet model
 *
 * This implementation follows the conventions for the naturally flavor
 * conserving 2HDM from [1106.0034](https://arxiv.org/abs/1106.0034).
 *
 */
class R2HDM : public TwoHDM {
public:
  //! short description
  static constexpr auto description = "CP-conserving 2HDM";
  //! neutral Higgs names, \f$ h, H, A\f$
  static constexpr std::array namesHzero{"Hl", "Hh", "A"};

  /**
   * @brief Input parametrization in terms of mixing angles.
   *
   * No ordering of \f$m_{H_a},m_{H_b}\f$ is assumed. Instead, the input mixing
   * matrix \f$R^\mathrm{in}\f$ in the basis \f$(H_a, H_b)\f$ is calculated
   * using the input mixing angle\f$\alpha^\mathrm{in}\f$. The mass ordered
   * states \f$h\f$ and \f$H\f$ and the corresponding mixing angle \f$alpha\f$
   * in the conventional \f$(H,h)\f$ basis are automatically obtained.
   */
  struct AngleInput {
    double mHa;   //!< \f$ m_{H_a} \f$, neutral CP-even Higgs mass in GeV
    double mHb;   //!< \f$ m_{H_b} \f$, neutral CP-even Higgs mass in GeV
    double mA;    //!< \f$ m_A \f$, neutral CP-odd Higgs mass in GeV
    double mHp;   //!< \f$ m_{H^\pm} \f$, charged Higgs mass in GeV
    double alpha; //!< input mixing angle \f$ \alpha^i \f$
    double tbeta; //!< \f$\tan\beta\f$, ratio of vacuum expectation values
    double m12sq; //!< \f$m_{12}^2\f$, in \f$\mathrm{GeV}^2\f$
    Yuk type;     //!< the Yukawa type
    double v;     //!< EW vev in GeV
  };

  /**
   * @brief Physical Input parametrization.
   *
   * A reparametrization using the gauge coupling \f$ c(H_y VV) \f$ instead of
   * the mixing angle. For the case \f$H_x\equiv h_{125}\f$ this parametrization
   * can be used to force \f$H_x\f$ to the alignment limit, by constraining the
   * gauge coupling of \f$H_y\f$ to be small.
   */
  struct PhysicalInput {
    double mHa;    //!< @copybrief AngleInput::mHa
    double mHb;    //!< @copybrief AngleInput::mHb
    double mA;     //!< @copybrief AngleInput::mA
    double mHp;    //!< @copybrief AngleInput::mHp
    double c_HbVV; //!< \f$ c(H_b VV) \f$, effective gauge coupling of \f$H_b\f$
    double tbeta;  //!< @copybrief AngleInput::tbeta
    double m12sq;  //!< @copybrief AngleInput::m12sq
    Yuk type;      //!< @copybrief AngleInput::type
    double v;      //!< @copybrief AngleInput::v
  };

  //! A R2HDM parameter point
  struct ParameterPoint {
    //! light neutral CP-even Higgs mass \f$m_h\f$ in GeV
    const double mHl;
    //! heavy neutral CP-even Higgs mass \f$m_H\f$ in GeV
    const double mHh;
    //! neutral CP-odd Higgs mass \f$m_A\f$ in GeV
    const double mA;
    //! charged Higgs mass \f$m_{H^\pm}\f$
    const double mHp;
    //! \f$\tan\beta\f$, ratio of vacuum expectation values
    const double tbeta;
    //! soft \f$Z_2\f$ breaking parameter \f$m_{12}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m12sq;
    //! mixing angle \f$\alpha\f$ of the CP-even sector
    const double alpha;
    //! quartic potential parameters \f$\lambda_{1,2,3,4,5}\f$
    const std::array<double, 5> L;
    //! doublet mass term \f$m_{11}^2\f$
    const double m11sq;
    //! doublet mass term \f$m_{22}^2\f$
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

    //! output names fo the parameters
    static constexpr std::array parameterNames{
        "mHl", "mHh", "mA", "mHp", "tbeta", "m12sq", "alpha",   "L1",
        "L2",  "L3",  "L4", "L5",  "m11sq", "m22sq", "yuktype", "v"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Model implementation for Constraints::BFB
   *
   * Uses the TwoHDM::BFB implementation.
   *
   * @param L the quartic parameters of the scalar potential
   * @return is the scalar potential at `p` bounded from below?
   */
  static bool BFB(const std::array<double, 5> &L);

  /**
   * @brief Model implementation for Constraints::Unitarity
   *
   * Uses the TwoHDM::MaxUnitarityEV implementation.
   *
   * @param L Quartic couplings
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(const std::array<double, 5> &L);

  /**
   * @brief Model implementation for Constraints::AbsoluteStability
   *
   * Implements eq (16) from [1303.5098](https://arxiv.org/abs/1303.5098).
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
   * p.data with the names given in AnyHdecay::Hdecay::r2hdmKeys.
   *
   * @param p the parameter point
   */
  static void RunHdecay(ParameterPoint &p);

  /**
   * @brief Calculates Higgs couplings.
   *
   * Calculates and stores the effective  couplings \f$c(H_iVV)\f$ (`c_HiVV`)
   * and \f$c(H_iAZ)\f$ (`c_HiAZ`) as well as the CP-even fermion couplings
   * \f$c^e(H_if\bar{f})\f$ `c_Hiuu_e', `c_Hidd_e`, `c_Hill_e` and the CP-odd
   * fermion couplings \f$c^o(Af\bar{f})\f$ `c_Auu_o`, `c_Add_o`, `c_All_o` for
   * `Hi` either `Hl` or `Hh`. Also calculates and stores the
   * \f$g_{H_iH^+H^-}\f$ couplings as `c_HiHpHm` in the convention of eq.
   * (B.12,B.13) from [1507.00933](https://arxiv.org/abs/1507.00933).
   *
   * @param p The parameter point.
   */
  static void CalcCouplings(ParameterPoint &p);

  /**
   * @brief Calculate and store some LHC production cross sections
   *
   * **Requires CalcCouplings to be called beforehand.**
   *
   * Stores the 13TeV LHC \f$gg\to H_i\f$ cross sections as `x_Hi_ggH` and the
   * \f$pp\to b\bar{b}H_i\f$ cross sections as `x_Hi_bbH`, for `Hi` in
   * namesHzero. Uses Tools::SushiTables.
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
   * \bar{t}bH^\pm\f$ cross section as `x_tHpm` (in pb).
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

  //! BSMPT model name for Constraints::EWPT
  static constexpr auto bsmptModelName = "r2hdm";

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
