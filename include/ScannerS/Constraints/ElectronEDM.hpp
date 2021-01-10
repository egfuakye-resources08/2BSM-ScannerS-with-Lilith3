#pragma once

#include "ScannerS/Constraints/Constraint.hpp"

namespace ScannerS::Constraints {

/**
 * @brief Constraint from the electron electric dipole moment in CP-violating
 * models.
 *
 * @tparam Model a model class with a `Model::CalcElectronEDM(p)->double`
 * function for a corresponding ParameterPoint `p`. The function returns the
 * value of the electron EDM in \f$ e\,\mathrm{cm} \f$.
 */
template <class Model>
class ElectronEDM : public Constraint<ElectronEDM, Model> {
public:
  static constexpr auto constraintId = "eEDM"; //!< unique constraint ID

  //! 90% c.l. limit by ACME, Nature 562 (2018) in \f$ e\,\mathrm{cm} \f$
  static constexpr double eEdmLimit = 1.1e-29;

  //! constructor that sets the severity
  explicit ElectronEDM(Severity severity)
      : Constraint<ElectronEDM, Model>{severity} {}

  /**
   * @brief Obtains the electron EDM bound.
   *
   * Stores the value of the electron EDM as `edm_e` (in \f$e\,\mathrm{cm}\f$).
   * The `Model::CalcElectronEDM` function may store additional values, eg
   * individual contributions to the EDM.
   *
   * @param p the parameter point
   * @return is the edm compatible with the EDM limit
   */
  bool Apply(typename Model::ParameterPoint &p) {
    auto result = Model::CalcElectronEDM(p);
    p.data.Store("edm_e", result);
    return std::abs(result) < eEdmLimit;
  }
};
} // namespace ScannerS::Constraints
