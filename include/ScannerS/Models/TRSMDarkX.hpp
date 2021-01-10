#pragma once

#include "ScannerS/DataMap.hpp"
#include "ScannerS/Models/TRSM.hpp"
#include <array>

namespace ScannerS {

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {
/**
 * @brief The dark X phase of the TRSM, *work in progress*
 */
class TRSMDarkX : public TRSM {
public:
  //! short description
  static constexpr auto description = "TRSM dark X phase";
  //! @cond

  //! neutral Higgs names, \f$ h, H, h_x \f$
  static constexpr std::array namesHzero{"Hl", "Hh", "HX"};

  /**
   * @brief Input parametrization in terms of mixing angles
   *
   * @todo description
   *
   */
  struct AngleInput {
    double mHa;   //!< \f$ M_a \f$, neutral Higgs mass in GeV
    double mHb;   //!< \f$ M_b \f$, neutral Higgs mass in GeV
    double mHX;   //!< \f$ M_X \f$, dark Higgs mass in GeV
    double a;     //!< \f$ a \f$, input mixing angle
    double v;     //!< EW vev in GeV
    double vs;    //!< \f$ v_S \f$ in GeV
    double lamX;  //!< \f$\lambda_X\f$, the \f$X^4\f$ self coupling
    double lamHX; //!< \f$\lambda_{\PhiX}\f$, the \f$\Phi^2 X^2\f$ coupling
    double lamSX; //!< \f$\lambda_{SX}\f$ the \f$ S^2 X^2 \f$ coupling
  };

  //! A dark-X-phase TRSM parameter point
  struct ParameterPoint {
    const double mHl;
    const double mHh;
    const double mHD;              //< DM mass
    const double alpha;            //< mixing angle
    const double v;                //< the EW vev
    const double vs;               //< the first real singlet vev
    const std::array<double, 6> L; //< quartic couplings: lambdaH, lambdaS,
                                   // lambdaX, lambdaHS, lambdaHX, lambdaSX
    const double muHsq;            //< doublet mass term
    const double muSsq;            //< first singlet mass term
    const double muXsq;            //< second singlet mass term
    DataMap data;                  //< place for all additional data

    static constexpr std::array parameterNames{
        "mH1", "mH2", "mHD", "alpha", "v",     "vs",    "LH",   "LS",
        "LX",  "LHS", "LHX", "LSX",   "muHsq", "muSsq", "muXsq"};

    explicit ParameterPoint(const AngleInput &in);
  };

  // ====== TRSMBroken_EWP ======
  /**
   * @brief Provides input for the general STU constraint.
   *
   * Provides the required input to the calculator of the oblique parameters.
   *
   * @param p the parameter point
   * @return Constraints::STUDetail::STUParameters
   */
  static Constraints::STUDetail::STUParameters
  STUInput(const ParameterPoint &p);

  //! @endcond
};

} // namespace Models
} // namespace ScannerS
