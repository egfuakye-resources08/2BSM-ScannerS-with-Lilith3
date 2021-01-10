#pragma once

#include "ScannerS/Constraints/Constraint.hpp"

namespace ScannerS::Constraints {

/**
 * @brief Constraint from the requirement of boundedness from below.
 *
 * Uses conditions for boundedness *in the strong sense*, ie strict positivity
 * of the scalar potential for large fields. The conditions have to be implemted
 * for each model.
 *
 * @tparam Model A model class with a function `Model::BFB(p.L)->bool` for a
 * corresponding ParameterPoint `p` with quartic couplings `p.L`. The function
 * should return true exactly if the scalar potential at `p` is bounded from
 * below.
 */
template <class Model> class BFB : public Constraint<BFB, Model> {
public:
  static constexpr auto constraintId = "BFB"; //!< unique constraint ID

  //! Constructor that sets the severity
  explicit BFB(Severity severity) : Constraint<BFB, Model>{severity} {}

  /**
   * @brief Obtains the BFB limit.
   *
   * Stores no output quantities.
   *
   * @param p the parameter point
   * @return is the scalar potential at `p` bounded from below?
   */
  bool Apply(const typename Model::ParameterPoint &p) {
    return Model::BFB(p.L);
  }
};

} // namespace ScannerS::Constraints
