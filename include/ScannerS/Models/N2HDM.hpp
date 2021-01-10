#pragma once

#include <Eigen/Core>
#include <array>

namespace EVADE {
namespace Models {
class N2HDM; // IWYU pragma: keep
}
} // namespace EVADE

namespace ScannerS {

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {

/**
 * @brief Base class for the different phases and variants of the N2HDM
 *
 * The implementation follows the conventions of
 * [1612.01309](https://arxiv.org/abs/1612.01309)
 */
class N2HDM {
public:
  //! number of neutral Higgs bosons
  static constexpr int nHzero = 4;
  //! number of charged Higgs bosons
  static constexpr int nHplus = 1;

  /**
   * @brief Model implementation for Constraints::BFB
   *
   * Implements Eq. (3.50) from [1612.01309](https://arxiv.org/abs/1612.01309).
   *
   * @param L the quartic parameters of the scalar potential
   * @return is the scalar potential at `p` bounded from below?
   */
  static bool BFB(const std::array<double, 8> &L);

  /**
   * @brief Model implementation for Constraints::Unitarity
   *
   * Implements Eq. (3.43-3.48) from
   * [1612.01309](https://arxiv.org/abs/1612.01309).
   *
   * @param L Quartic couplings
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(const std::array<double, 8> &L);

  /**
   * @brief Provides input for the general STU calculation
   *
   * @param mA \f$ m_A \f$
   * @param mHi the three CP-even Higgs masses \f$ m_{H_{1,2,3}}\f$
   * @param mHp charged Higgs mass \f$m_{H^\pm}\f$
   * @param tbeta \f$\tan\beta\f$
   * @param R 3x3 CP-even mixing matrix
   * @return Constraints::STUDetail::STUParameters input parameters for the STU
   * calculation
   */
  static Constraints::STUDetail::STUParameters
  STUInput(double mA, const std::array<double, nHzero - 1> &mHi, double mHp,
           double tbeta, const Eigen::Matrix3d &R);

  //! EVADE model for Constraints::VacStab
  using ModelEVADE = EVADE::Models::N2HDM;
};

} // namespace Models
} // namespace ScannerS
