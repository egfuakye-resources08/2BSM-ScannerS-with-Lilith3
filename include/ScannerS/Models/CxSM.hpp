#pragma once

#include <array>

namespace ScannerS::Models {
/**
 * @brief Base class for the different phases of the CxSM
 *
 * The implementation follows the conventions of
 * [1512.05355](https://arxiv.org/abs/1512.05355)
 */
class CxSM {
public:
  //! number of neutral Higgs bosons
  static constexpr int nHzero = 3;
  //! number of charged Higgs bosons
  static constexpr int nHplus = 0;
  //! charged Higgs names, empty
  static constexpr std::array<const char *, 0> namesHplus{};

  /**
   * @brief Model implementation for Constraints::BFB
   *
   * Implements Eq. (3.11) from [1301.2599](https://arxiv.org/abs/1301.2599).
   *
   * @param L Quartic couplings
   * @return bounded?
   */
  static bool BFB(const std::array<double, 3> &L);

  /**
   * @brief Model implementation for Constraints::Unitarity
   * @param L quartic couplings
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(const std::array<double, 3> &L);

  /**
   * @brief Model implementation for Constraints::STU
   * @return `true`, the oblique parameter approximation is always valid for
   * singlet extensions
   */
  template <class T> static bool EWPValid(const T &) { return true; }
};
} // namespace ScannerS::Models
