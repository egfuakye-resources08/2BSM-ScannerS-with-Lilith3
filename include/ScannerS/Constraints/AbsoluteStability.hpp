#pragma once

#include "ScannerS/Constraints/Constraint.hpp"

namespace ScannerS::Constraints {

/**
 * @brief Constraint to require absolute stability of the EW vacuum.
 *
 * This is meant to implement conditions like the abslute stability discriminant
 * in the R2HDM from [1303.5098](https://arxiv.org/abs/1303.5098).
 *
 * @tparam Model a model class that provides a function
 * `Model::AbsoluteStability(p)->bool` for a corresponding parameter point `p`.
 * The functions should return whether the EW vacuum at `p` is absolutely
 * stable.
 */
template <class Model>
class AbsoluteStability : public Constraint<AbsoluteStability, Model> {
public:
  static constexpr auto constraintId = "AbsStab"; //!< unique constraint ID

  //! Constructor that sets the severity
  explicit AbsoluteStability(Severity severity)
      : Constraint<AbsoluteStability, Model>{severity} {}

  /**
   * @brief Obtains the absolute stability bound
   *
   * Stores no output quantities.
   *
   * @param p the parameter point
   * @return is the EW vacuum at `p` absolutely stable
   */
  bool Apply(const typename Model::ParameterPoint &p) {
    return Model::AbsoluteStability(p);
  }
};

} // namespace ScannerS::Constraints
