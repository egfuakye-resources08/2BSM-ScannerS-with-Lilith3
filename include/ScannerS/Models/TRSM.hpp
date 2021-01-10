#pragma once

#include "ScannerS/DataMap.hpp"
#include <Eigen/Core>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace EVADE {
namespace Models {
class TRSM; // IWYU pragma: keep
}
} // namespace EVADE

namespace ScannerS {

namespace Constraints::STUDetail {
struct STUParameters;
}

namespace Models {

/**
 * @brief Base class for the different phases of the Two-Real-Singlet Model
 *
 * The implementation follows the conventions of
 * [1908.08554](https://arxiv.org/abs/1908.08554)
 */
class TRSM {
public:
  //! number of neutral Higgs bosons
  static constexpr int nHzero = 3;
  //! number of charged Higgs bosons
  static constexpr int nHplus = 0;
  //! charged Higgs names, empty
  static constexpr std::array<const char *, nHplus> namesHplus{};

  /**
   * @brief Model implementation for Constraints::BFB
   *
   * Implements Eq. (62) from [1603.02680](https://arxiv.org/abs/1603.02680).
   *
   * @param L the quartic parameters of the scalar potential
   * @return is the scalar potential at `p` bounded from below?
   */
  static bool BFB(const std::array<double, 6> &L);

  /**
   * @brief Model implementation for Constraints::Unitarity
   *
   * Implements eqs. (30-32) from [1908.08554](https://arxiv.org/abs/1908.08554)
   *
   * @param L Quartic couplings
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(const std::array<double, 6> &L);

  /**
   * @brief Provides input for the general STU calculation
   *
   * @param mHi the three TRSM Higgs masses
   * @param Ri0 doublet mixing elements of the three scalars
   * @return Constraints::STUDetail::STUParameters input parameters for the STU
   * calculation
   */
  static Constraints::STUDetail::STUParameters
  STUInput(const std::array<double, nHzero> &mHi,
           const std::array<double, nHzero> &Ri0);

  /**
   * @brief Model implementation for Constraints::STU
   * @return `true`, the oblique parameter approximation is always valid for
   * singlet extensions
   */
  template <class T> static bool EWPValid(const T &) { return true; }

  /**
   * @brief Calculate all TRSM triple Higgs couplings \f$\tilde\lambda_{abc}\f$.
   *
   * Phase-independent implementation in terms of lagrangian parmeters, a 3x3
   * mixing matrix (that may be block-diagonal in the dark phases) and vevs.
   * Uses the conventions of eqs. (19-21) from
   * [1908.08554](https://arxiv.org/abs/1908.08554)
   *
   * @param L quartic potential parameters
   * @param R mixing matrix
   * @param v doublet vev
   * @param vs singlet vev
   * @param vx singlet vev
   * @return a map of the triple Higgs couplings with keys of the form
   * `c_HaHbHc`
   */
  static DataMap::Map TripleHCoups(const std::array<double, 6> &L,
                                   const Eigen::Matrix3d &R, double v,
                                   double vs, double vx);

  /**
   * @brief Model implementation for Constraints::VacStab
   * @param p the parameter point
   * @return std::vector<double> vector of EVADE input parameters
   */
  template <class ParamP>
  static std::vector<double> ParamsEVADE(const ParamP &p) {
    return {p.muHsq, p.muSsq, p.muXsq, p.L[0], p.L[1], p.L[2],
            p.L[3],  p.L[4],  p.L[5],  p.v,    p.vs,   p.vx};
  }

  //! EVADE model for Constraints::VacStab
  using ModelEVADE = EVADE::Models::TRSM;
};

} // namespace Models
} // namespace ScannerS
