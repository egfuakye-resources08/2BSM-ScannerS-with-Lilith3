#pragma once

#include "ScannerS/Constants.hpp"
#include "ScannerS/Constraints/Constraint.hpp"

namespace ScannerS::Constraints {

/**
 * @brief Constraint from tree-level perturbative unitarity.
 *
 * Sets an upper limit on the absolute value of the largest eigenvalue of the
 * \f$ 2\to2 \f$ scatter matrix \f$ \mathcal{M}_{2\to2} \f$. The eigenvalues
 * have to be implemted for each model.
 *
 * @tparam Model A model class with a function
 * `Model::MaxUnitarityEV(p.L)->double` for a corresponding ParameterPoint `p`
 * with quartic couplings `p.L`. The function should return the largest absolute
 * value of the eigenvalues of the scattering matrix.
 */
template <class Model> class Unitarity : public Constraint<Unitarity, Model> {
public:
  static constexpr auto constraintId = "Uni"; //!< unique constraint ID

  //! Constructor that sets the severity and the upper limit on maxEV
  explicit Unitarity(Severity severity,
                     double unitarityLimit = 8 * Constants::pi)
      : Constraint<Unitarity, Model>{severity}, unitarityLimit_{
                                                    unitarityLimit} {}

  /**
   * @brief Obtains the unitarity limit.
   *
   * Stores \f$ \max(|\mathcal{M}^i_{2\to2}|) \f$ as `maxEV`.
   *
   * @param p the parameter point
   * @return is `p` unitary for the given `unitarityLimit` on `maxEV`
   */
  bool Apply(typename Model::ParameterPoint &p) {
    double maxEV = Model::MaxUnitarityEV(p.L);
    p.data.Store("maxEV", maxEV);
    return maxEV < unitarityLimit_;
  }

private:
  const double unitarityLimit_;
};

} // namespace ScannerS::Constraints
