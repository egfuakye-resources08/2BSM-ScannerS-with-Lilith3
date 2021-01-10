#pragma once

#include "ScannerS/Constraints/Constraint.hpp"
#include <stdexcept>

namespace ScannerS::Constraints {

//! functions related to the BPhysics constraint
namespace BPhysicsDetail {

/**
 * @brief 2HDM T1 \f$ B_d \to \mu \mu \f$ constraint.
 * From 1803.01853 Fig 9. Lower bound on tbeta.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool T1Bdmumu(double tbeta, double mHp);

/**
 * @brief 2HDM T2 \f$ B_s \to \mu \mu \f$ constraint.
 * From 1803.01853 Fig 9. Upper bound on tbeta.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool T2Bsmumu(double tbeta, double mHp);

/**
 * @brief 2HDM T2 \f$ B \to X_s \gamma \f$ constraint.
 * From 1803.01853 Fig 9. Lower bound on mHp.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool T2Bsgam(double tbeta, double mHp);

/**
 * @brief 2HDM lepton specific \f$ B_d \to \mu \mu \f$ constraint.
 * From 1803.01853 Fig 9. Lower bound on tbeta.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool LSBdmumu(double tbeta, double mHp);

/**
 * @brief 2HDM flipped \f$ B_d \to \mu \mu \f$ constraint.
 * From 1803.01853 Fig 9. Lower bound on tbeta.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool FBdmumu(double tbeta, double mHp);

/**
 * @brief 2HDM flipped \f$ B \to X_s \gamma \f$ constraint.
 * From 1803.01853 Fig 9. Lower bound on mHp.
 * @relatedalso ScannerS::Constraints::BPhysics
 */
bool FBsgam(double tbeta, double mHp);
} // namespace BPhysicsDetail

/**
 * @brief Constraint from B-physics flavor observables.
 *
 * Only applicable to 2HDM-like models with a Yukawa \f$ \mathbb{Z}_2 \f$
 * symmetry and a single charged Higgs boson. Uses the fit results from Fig. 9
 * of [1803.01853](https://arxiv.org/abs/1803.01853).
 *
 * @tparam Model a model class that has a `Model::Yuk` enumeration member type
 * equivalent to 2HDM::Yuk. For a corresponding ParameterPoint `p`, `p.type` has
 * to be a member of this enumeration, `p.mHp` has to be the charged Higgs mass,
 * and `p.tbeta` has to correspond to the 2HDM parameter \f$ \tan\beta \f$ as
 * far as the charged Higgs couplings to fermions are concerned.
 */
template <class Model> class BPhysics : public Constraint<BPhysics, Model> {
public:
  static constexpr auto constraintId = "BPhys"; //!< unique constraint ID

  //! Constructor that sets the severity
  explicit BPhysics(Severity severity)
      : Constraint<BPhysics, Model>{severity} {}

  /**
   * @brief Obtains the bound from B-physics flavor observables.
   *
   * Stores no output quantities.
   *
   * @param p the parameter point
   * @return are the B-physics constraints fulfilled
   */
  bool Apply(const typename Model::ParameterPoint &p) {
    switch (p.type) {
    case Model::Yuk::typeI:
      return BPhysicsDetail::T1Bdmumu(p.tbeta, p.mHp);
    case Model::Yuk::typeII:
      return BPhysicsDetail::T2Bsmumu(p.tbeta, p.mHp) &&
             BPhysicsDetail::T2Bsgam(p.tbeta, p.mHp);
    case Model::Yuk::leptonSpecific:
      return BPhysicsDetail::LSBdmumu(p.tbeta, p.mHp);
    case Model::Yuk::flipped:
      return BPhysicsDetail::FBdmumu(p.tbeta, p.mHp) &&
             BPhysicsDetail::FBsgam(p.tbeta, p.mHp);
    }
    throw(std::runtime_error("Invalid Yukawa type."));
  }
};
} // namespace ScannerS::Constraints
