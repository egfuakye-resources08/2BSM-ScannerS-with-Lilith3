#pragma once

#include "ScannerS/DataMap.hpp"
#include <BSMPT/models/ClassPotentialOrigin.h> // IWYU pragma: keep
#include <BSMPT/models/IncludeAllModels.h>
#include <ScannerS/Constraints/Constraint.hpp>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ScannerS::Constraints {

//! additional functions for the EWPT constraint
namespace EWPTDetail {

//! initialize the EWPT model and return the counterterms
//! @relatedalso ScannerS::Constraints::EWPT
DataMap::Map InitModel(std::shared_ptr<BSMPT::Class_Potential_Origin> model,
                       const std::vector<double> &input);

//! obtain the triple Higgs couplings
//! @relatedalso ScannerS::Constraints::EWPT
DataMap::Map
TripleHiggsCoups(std::shared_ptr<BSMPT::Class_Potential_Origin> model);

//! find the NLO T=0 VEV and check if it matches the EW vacuum
//! @relatedalso ScannerS::Constraints::EWPT
std::pair<DataMap::Map, bool>
NloVev(std::shared_ptr<BSMPT::Class_Potential_Origin> model);

//! Calculate the EW phase transition and return its properties
//! @relatedalso ScannerS::Constraints::EWPT
DataMap::Map
PhaseTransition(std::shared_ptr<BSMPT::Class_Potential_Origin> model);
} // namespace EWPTDetail

/**
 * @brief Constraint from the requirement of a first order EW phase transition.
 *
 * This is not a true constraint, but rather a requirement if gravitational wave
 * signatures or EW baryogenesis should be studied.
 *
 * Uses the code [BSMPT](https://github.com/phbasler/BSMPT) to perform the
 * calculation. **Relies on BSMPT-2.0 or newer, which is not yet public.***
 *
 * @tparam Model a model class that provides a
 * `Model::BsmptInput(p)->std::vector<double>` function that returns the input
 * parameters for BSMPT for a given ParameterPoint `p`. The model must also
 * provide a `Model::bsmptModelName` member that matches the corresponding name
 * in `BSMPT/models/IncludeAllModels.h`.
 */
template <class Model> class EWPT : public Constraint<EWPT, Model> {
public:
  static constexpr auto constraintId = "EWPT"; //!< unique constraint ID

  /**
   * @brief Constructor
   *
   * @param severity the severity
   * @param minimumPtStrength minimum required strength of the first order phase
   * transition
   */
  explicit EWPT(Severity severity, double minimumPtStrength = 0.)
      : Constraint<EWPT, Model>{severity}, minimumPtStrength_{
                                               minimumPtStrength} {}

  /**
   * @brief Obtains the bound from requiring a first order EWPT
   *
   * Stores *lots* of values, including all counterterms, all triple Higgs
   * couplings, the NLO T=0 vacuum, and the field configuration at the critical
   * PT temperature. The names are model specific and documented in BSMPT.
   * Regardless of the model, this always stores the EWPT vev as `EWPT_omega_c`,
   * the critical temperature as `EWPT_T_c` and a status flag in `BSMPT_ok`.
   *
   * @param p the parameter point
   * @return is the EWPT first order, and the T=0 vev absolutely stable at NLO
   */
  bool Apply(typename Model::ParameterPoint &p) {
    auto input = Model::BsmptInput(p);
    p.data.Merge(EWPTDetail::InitModel(bsmptModel, input));
    p.data.Merge(EWPTDetail::TripleHiggsCoups(bsmptModel));
    auto [nloVev, nloStable] = EWPTDetail::NloVev(bsmptModel);
    p.data.Merge(std::move(nloVev));
    p.data.Merge(EWPTDetail::PhaseTransition(bsmptModel));

    return nloStable && p.data["BSMPT_ok"] > 0 &&
           p.data["EWPT_omega_c"] > minimumPtStrength_ * p.data["EWPT_T_c"];
  }

private:
  const std::shared_ptr<BSMPT::Class_Potential_Origin> bsmptModel =
      BSMPT::ModelID::FChoose(BSMPT::ModelID::getModel(Model::bsmptModelName));
  double minimumPtStrength_;
};
} // namespace ScannerS::Constraints
