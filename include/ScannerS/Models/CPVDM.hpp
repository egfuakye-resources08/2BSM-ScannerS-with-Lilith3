#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/N2HDM.hpp"
#include <Eigen/Core>
#include <array>
#include <complex>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace EVADE {
namespace Models {
class CDN2HDM; // IWYU pragma: keep
}
} // namespace EVADE

namespace ScannerS {
namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Interfaces {
namespace HiggsBoundsSignals {
template <int nHzero, int nHplus> struct HBInputEffC;
template <size_t nHzero, size_t nHplus> class HiggsBoundsSignals;
} // namespace HiggsBoundsSignals
} // namespace Interfaces

namespace Models {
/**
 * @brief The minimal model of CP-violating scalar dark matter
 *
 * The implementation follows [1807.10322](https://arxiv.org/abs/1807.10322).
 *
 */
class CPVDM : public N2HDM {
public:
  //! short description
  static constexpr auto description = "minimal CP-violating dark matter model";
  //! neutral Higgs names, \f$ h, h_1, h_2, h_3\f$
  static constexpr std::array namesHzero{"Hsm", "H1", "H2", "H3"};
  //! charged Higgs names, \f$ H^\pm \f$
  static constexpr std::array namesHplus{"Hp"};

  /**
   * @brief Input parametrization in terms of mixing angles
   *
   * The third dark neutral Higgs mass \f$m_{h_c}\f$ is calculated from
   * \f$m_{h_a}, m_{h_b}\f$ and the initial mixing matrix \f$ R^\mathrm{in} \f$
   * in the basis \f$h_a,h_b,h_c\f$ (= Utilities::MixMat3d() with the input
   * mixing angles \f$ \alpha^\mathrm{in}_{1,2,3}\f$). No ordering of \f$
   * m_{h_a}, m_{h_b}, m_{h_c}\f$ is required. The mass-ordered mixing matrix
   * \f$ R \f$ and the corresponding mass-ordered eigenstates \f$ h_{1,2,3} \f$
   * and mixing angles \f$\alpha_{1,2,3}\f$ are automatically obtained.
   *
   */
  struct AngleInput {
    double mHsm;  //!< \f$ m_h \f$, visible Higgs mass in GeV
    double mHa;   //!< \f$ m_{h_a} \f$, dark neutral Higgs mass in GeV
    double mHb;   //!< \f$ m_{h_b} \f$, dark neutral Higgs mass in GeV
    double mHp;   //!< \f$m_{H^\pm}\f$, dark charged Higgs mass in GeV
    double a1;    //!< first mixing angle \f$\alpha^\mathrm{in}_1\f$
    double a2;    //!< second mixing angle \f$\alpha^\mathrm{in}_2\f$
    double a3;    //!< third mixing angle \f$\alpha^\mathrm{in}_3\f$
    double L2;    //!< \f$\lambda_2\f$
    double L6;    //!< \f$\lambda_6\f$
    double L8;    //!< \f$\lambda_8\f$
    double m22sq; //!< \f$m_{22}^2\f$, in \f$\mathrm{GeV}^2\f$
    double mssq;  //!< \f$ m_S^2 \f$, in \f$\mathrm{GeV}^2\f$
    double v;     //!< EW vev in GeV
  };

  //! Parameter point of the CPVDM
  struct ParameterPoint {
    //! visible Higgs mass \f$ m_h \f$ in GeV
    const double mHsm;
    //! mass-ordered dark neutral Higgs masses \f$ m_{h_{1,2,3}}\f$ in GeV
    const std::array<double, 3> mHi;
    //! dark charged Higgs mass \f$ m_{H^\pm} \f$ in GeV
    const double mHp;
    //! mass-ordered dark neutral Higgs mixing matrix \f$R\f$
    const Eigen::Matrix3d R;
    //! quartic potential parameters \f$\lambda_{1,2,3,4,5,6,7,8}\f$
    const std::array<double, 8> L;
    //! trilinear parameter \f$A\f$
    const std::complex<double> A;
    //! doublet mass term \f$m_{11}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m11sq;
    //! doublet mass term \f$m_{22}^2\f$ in \f$\mathrm{GeV}^2\f$
    const double m22sq;
    //! singlet mass term \f$m_S^2\f$ in \f$\mathrm{GeV}^2\f$
    const double mssq;
    //! EW vev in GeV
    const double v;
    //! dark sector mixing angles \f$\alpha_{1,2,3}\f$ parametrizing \f$R\f$
    const std::array<double, 3> alpha;
    //! place for additional data
    DataMap data;

    //! Construct a ParameterPoint from AngleInput
    explicit ParameterPoint(const AngleInput &in);

    //! output names for the parameters
    static constexpr std::array parameterNames{
        "mHsm",  "mH1",   "mH2",  "mH3", "mHp", "R11", "R12", "R13",
        "R21",   "R22",   "R23",  "R31", "R32", "R33", "L1",  "L2",
        "L3",    "L4",    "L5",   "L6",  "L7",  "L8",  "Tr",  "Ti",
        "m11sq", "m22sq", "mssq", "v",   "a1",  "a2",  "a3"};

    //! serialize the parameter and data values for output
    std::string ToString() const;
  };

  /**
   * @brief Checks if the generation in the ParameterPoint constructor was
   * successful, **ALWAYS CALL THIS**.
   *
   * The construction of a ParameterPoint in the CPVDM is not guaranteed to work
   * for all possible inputs. Since the third neutral Higgs mass is calculated
   * from the other two and the mixing matrix, it is not guaranteed to be
   * physical. If the input would result in a tachyonic third Higgs mass,
   * \f$m^2_{h_z}<0\f$, it is set to `-1`. Since `p.mHi` is mass-ordered, this
   * means that `p` is invalid if \f$m_{h_1}<0\f$.
   *
   * @param p the parameter point
   * @return is `p` valid
   */
  inline static bool Valid(const ParameterPoint &p) { return p.mHi[0] > 0; }

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
   * @brief Calculates Higgs couplings.
   *
   * Calculates and stores all triple Higgs couplings between \f$h\f$ and
   * \f$h_{1,2,3}\f$ as `c_HsmHiHj` for `i,j` in `1,2,3` as well as the coupling
   * between \f$h H^\pm H^\mp \f$ as `c_HsmHpHm`.
   *
   * @param p the parameter point
   */
  static void CalcCouplings(ParameterPoint &p);

  /**
   * @brief Model implementation for Constraints::Higgs
   *
   * **Requires CalcCouplings to be called beforehand.**
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
   * @brief Model implementation for Constraints::VacStab
   * @param p the parameter point
   * @return std::vector<double> vector of EVADE input parameters
   */
  static std::vector<double> ParamsEVADE(const ParameterPoint &p);
  //! EVADE model for Constraints::VacStab
  using ModelEVADE = EVADE::Models::CDN2HDM;

  /**
   * @brief Model implementation for Constraints::DM
   * @param p the parameter point
   * @return std::map<std::string, double> of MicrOMEGAs input parameters
   */
  static std::map<std::string, double> MOInput(const ParameterPoint &p);
  //! Constraints::DM MicrOMEGAs model name
  static constexpr auto micromegasModelName = "CPVDM";
};

} // namespace Models
} // namespace ScannerS
