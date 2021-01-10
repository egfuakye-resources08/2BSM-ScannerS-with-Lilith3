#pragma once

#include <array>

namespace ScannerS::Models {

//! Base class for two Higgs doublet models
class TwoHDM {
public:
  //! number of neutral Higgs bosons
  static constexpr int nHzero = 3;
  //! number of charged Higgs bosons
  static constexpr int nHplus = 1;
  //! charged Higgs names, \f$ H^\pm \f$
  static constexpr std::array namesHplus{"Hp"};

  //! \f$Z_2\f$ Yukawa types
  enum class Yuk { typeI = 1, typeII = 2, leptonSpecific = 3, flipped = 4 };

  /**
   * @brief Checks if the scalar potential is bounded.
   *
   * Implements eq (176) from [1106.0034](https://arxiv.org/abs/1106.0034).
   *
   * @param L1 \f$ \lambda_1 \f$
   * @param L2 \f$ \lambda_2 \f$
   * @param L3 \f$ \lambda_3 \f$
   * @param L4 \f$ \lambda_4 \f$
   * @param abs_L5 \f$ |\lambda_5| \f$
   * @return is the potential bounded?
   */
  static bool BFB(double L1, double L2, double L3, double L4, double abs_L5);

  /**
   * @brief Calculates the largest eigenvalue of the tree-level \f$2\to2\f$
   * scattering matrix.
   *
   * Implements eq (372) from [1106.0034](https://arxiv.org/abs/1106.0034).
   *
   * @param L1 \f$ \lambda_1 \f$
   * @param L2 \f$ \lambda_2 \f$
   * @param L3 \f$ \lambda_3 \f$
   * @param L4 \f$ \lambda_4 \f$
   * @param abs_L5 \f$ |\lambda_5| \f$
   * @return maxEV, absolute value of largest eigenvalue
   */
  static double MaxUnitarityEV(double L1, double L2, double L3, double L4,
                               double abs_L5);

  //! Effective charged Higgs couplings to quarks
  struct HpCoups {
    double rhot; //!< effective coupling to top quarks
    double rhob; //!< effective coupling to bottom quarks
  };

  //! Charged Higgs couplings from \f$\tan\beta\f$ and the Yukawa type
  static HpCoups TwoHDMHpCoups(double tbeta, TwoHDM::Yuk type);
};

} // namespace ScannerS::Models
