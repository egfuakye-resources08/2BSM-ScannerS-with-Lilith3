#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/N2HDM.hpp"
#include "ScannerS/Models/TwoHDM.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ScannerS {
namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Tools {
class SushiTables;
} // namespace Tools

namespace Interfaces {
namespace HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInput;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace HiggsBoundsSignals
} // namespace Interfaces

namespace Models {
/**
 * @brief The broken phase of the N2HDM
 *
 * The implementation follows
 * [1612.01309](https://arxiv.org/abs/1612.01309).
 */
class N2HDMBroken : public N2HDM {
public:
  //! short description
  static constexpr auto description = "N2HDM broken phase";
  //! names of the neutral Higgs bosons, \f$ A, H_1, H_2, H_3 \f$
  static constexpr std::array namesHzero{"A", "H1", "H2", "H3"};
  //! charged Higgs names, \f$ H^\pm \f$
  static constexpr std::array namesHplus{"Hp"};
  //! Yukawa types
  using Yuk = TwoHDM::Yuk;

  /**
   * @brief Input parametrization in terms of mixing angles.
   *
   * No mass ordering between \f$m_{H_a},m_{H_b},m_{H_c}\f$ is required. The
   * initial mixing matrix \f$ R^\mathrm{in} \f$ corresponds to the ordering
   * \f$(H_a, H_b, H_c)\f$ and is calculated from the initial mixing angles
   * \f$\alpha^\mathrm{in}_{1,2,3}\f$ (using Utilities::MixMat3d()). The
   * mass-ordered mixing matrix \f$ R \f$ and the corresponding mass-ordered
   * eigenstates \f$ H_{1,2,3} \f$ and mixing angles \f$\alpha_{1,2,3}\f$ are
   * automatically obtained.
   */
  struct AngleInput {
    double mHa;   //!< \f$ m_{H_a} \f$, neutral CP-even Higgs mass in GeV
    double mHb;   //!< \f$ m_{H_b} \f$, neutral CP-even Higgs mass in GeV
    double mHc;   //!< \f$ m_{H_c} \f$, neutral CP-even Higgs mass in GeV
    double mA;    //!< \f$ m_A \f$, neutral CP-odd Higgs mass in GeV
    double mHp;   //!< \f$ m_{H^\pm} \f$, charged Higgs mass in GeV
    double tbeta; //!< \f$\tan\beta\f$, ratio of vacuum expectation values
    double a1;    //!< first mixing angle \f$\alpha^\mathrm{in}_1\f$
    double a2;    //!< second mixing angle \f$\alpha^\mathrm{in}_2\f$
    double a3;    //!< third mixing angle \f$\alpha^\mathrm{in}_3\f$
    double m12sq; //!< \f$ m_{12}^2\f$ in \f$\mathrm{GeV}^2\f$
    Yuk type;     //!< the Yukawa type
    double vs;    //!< singlet vev \f$ v_s \f$ in GeV
    double v;     //!< EW vev in GeV
  };

  /**
   * @brief Physical input parametrization
   *
   * A reparametrization using couplings instead of mixing angles to allow a
   * more straightforward selection of scan ranges compared to AngleInput when
   * considering \f$H_a\equiv h_{125}\f$.
   *
   * The AngleInput parameters \f$\alpha^\mathrm{in}_{1,2,3}\f$ are expressed
   * through \f$c^2(H_aVV)\f$, \f$c^2(H_at\bar{t})\f$, \f$R^\mathrm{in}_{b3}\f$
   * and \f$\mathrm{sign}(R^\mathrm{in}_{a1}\cdot R^\mathrm{in}_{a3})\f$. Two
   * signs are needed as input to compensate for the *squared* input for the
   * effective couplings. The first one is
   * \f$\mathrm{sign}(R^\mathrm{in}_{a1}\cdot R^\mathrm{in}_{a3})\f$ and the
   * second one is fixed through the *physical* assumption, that
   * \f$c(H_aVV)\times c(H_at\bar{t})>0\f$, which is a known consequence of
   * Higgs measurements for \f$H_a = h_{125}\f$.
   */
  struct PhysicalInput {
    double mHa;   //!< @copybrief AngleInput::mHa
    double mHb;   //!< @copybrief AngleInput::mHb
    double mHc;   //!< @copybrief AngleInput::mHc
    double mA;    //!< @copybrief AngleInput::mA
    double mHp;   //!< @copybrief AngleInput::mHp
    double tbeta; //!< @copybrief AngleInput::tbeta
    //! \f$c^2(H_aVV)\f$, squared effective coupling of \f$ H_a\f$ to W/Z bosons
    double c_HaVV_sq;
    //! \f$c^2(H_at\bar{t})\f$, squared effective coupling of \f$ H_a\f$ to
    //! top-quarks
    double c_Hatt_sq;
    //! \f$\mathrm{sign}(R^\mathrm{in}_{a3})\f$, sign of the mixing matrix
    //! element \f$R^\mathrm{in}_{a3}\f$
    int sign_Ra3;
    //! \f$R_{b3}\f$ mixing element between \f$H_b\f$ and the singlet field
    double Rb3;
    double m12sq; //!< @copybrief AngleInput::m12sq
    Yuk type;     //!< @copybrief AngleInput::type
    double vs;    //!< @copybrief AngleInput::vs
    double v;     //!< @copybrief AngleInput::v
  };

  //! Parameter point of the broken phase N2HDM
  struct ParameterPoint {
    //! mass-ordered neutral CP-even Higgs masses \f$m_{H_{1,2,3}}\f$ in GeV
    const std::array<double, nHzero - 1> mHi;
    //! CP-odd Higgs mass \f$m_A\f$ in GeV
    const double mA;
    //! charged Higgs mass \f$m_{H^\pm}\f$ in GeV
    const double mHp;
    //! \f$\tan\beta\f$, ratio of doublet vacuum expectation values
    const double tbeta;
    //! mass-ordered neutral mixing matrix \f$R\f$
    const Eigen::Matrix3d R;
    //! the corresponding mixing angles \f$\alpha_{1,2,3}\f$
    const std::array<double, 3> alpha;
    //! Yukawa type
    const Yuk type;
    //! singlet vev \f$ v_s \f$ in GeV
    const double vs;
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
    //! Construct a ParameterPoint from PhysicalInput
    explicit ParameterPoint(const PhysicalInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mH1",     "mH2", "mH3", "mA",    "mHp",   "tbeta", "R11", "R12", "R13",
        "R21",     "R22", "R23", "R31",   "R32",   "R33",   "a1",  "a2",  "a3",
        "yuktype", "vs",  "v",   "m12sq", "L1",    "L2",    "L3",  "L4",  "L5",
        "L6",      "L7",  "L8",  "m11sq", "m22sq", "mssq"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Checks if the generation in the ParameterPoint constructor was
   * successful, **ALWAYS CALL THIS**.
   *
   * The construction of a ParameterPoint from PhysicalInput in the N2HDM is not
   * guaranteed to work for all possible inputs. There are combinations of input
   * couplings for which no solution for \f$a_1,a_2,a_3\f$ exist. In this case,
   * the mixing matrix is set to a constant value \f$> 1\f$ such that `p` is
   * invalid if \f$R_{11}>1\f$.
   *
   * @param p the parameter point
   * @return is `p` valid
   */
  static inline bool Valid(const ParameterPoint &p) { return p.R(0, 0) <= 1; }

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
   * p.data with the names given in AnyHdecay::Hdecay::n2hdmBrokenKeys.
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
   * `c_Auu_o`, `c_Add_o`, `c_All_o`) couplings for `Hi` in `{"H1","H2","H3"}`.
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
   * @return std::vector<double> vector of EVADE input parameters
   */
  static std::vector<double> ParamsEVADE(const ParameterPoint &p);

  //! BSMPT model name for Constraints::EWPT
  static constexpr auto bsmptModelName = "n2hdm";

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
