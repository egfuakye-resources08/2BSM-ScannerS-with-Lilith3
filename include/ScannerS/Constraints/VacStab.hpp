#pragma once

#include "EVADE/EVADE.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace ScannerS::Constraints {

/**
 * @brief Constraint from metastability of the electroweak vacuum.
 *
 * Uses the EVADE library to determine all minima of the tree-level scalar
 * potential and calculate the tunnelling times to all possible deeper
 * stationary points.
 *
 * @tparam Model a model class that specifies its corresponding EVADE model
 * class through the `Model::ModelEVADE` type and provides a
 * `Model::ParamsEVADE(p)->std::vector<double>` function that returns the input
 * parameters expected by the EVADE model.
 */
template <class Model> class VacStab : public Constraint<VacStab, Model> {
public:
  static constexpr auto constraintId = "vacstab"; //!< unique constraint ID

  /**
   * @brief Constructor
   *
   * @param severity the severity
   * @param fieldSets the combinations of fields to consider for simultaneous
   * solutions to the stationarity conditions
   */
  VacStab(Severity severity, std::vector<std::vector<std::string>> fieldSets)
      : Constraint<VacStab, Model>{severity}, fieldSets_{std::move(fieldSets)},
        statConds_{EVADE::Solver::Hom4ps2(".")} {}

  /**
   * @brief Obtains the metastability bound
   *
   * Stores information on both the most dangerous minimum (prefixed `MDM_`) and
   * the flobal minimum (`GM_`). This information includes the value of the
   * scalar potential (`V`), the bounce action for tunnelling from the EW vacuum
   * (`B`) and the values of all field components at the respective stationary
   * point (see the FieldNames() function of the corresponding EVADE model
   * class).
   *
   * @param p the parameter point
   * @return is the EW vacuum at `p` at least metastable and long-lived
   */
  bool Apply(typename Model::ParameterPoint &p) {
    auto pars = Model::ParamsEVADE(p);
    auto statpoints = EVADE::SolveForFieldSets(statConds_, pars, fieldSets_);
    auto result =
        EVADE::CalculateTunnelling<EModel>(statpoints, pars); // pass results
    const auto fieldNames = EModel::FieldNames();

    auto storeFSP = [&fieldNames,
                     &p](const EVADE::Fieldspace::TunnellingDir &dir,
                         const std::string &key) {
      p.data.Store(key + "_B", dir.B());
      const auto fields = dir.Target().Fields();
      for (size_t i = 0; i != fieldNames.size(); ++i) {
        p.data.Store(key + "_" + fieldNames[i], fields[i]);
      }
      p.data.Store(key + "_V", dir.Target().V());
    };

    // store MDM
    const auto mdm = result[0];
    storeFSP(mdm, "MDM");

    // store global minimum
    auto lessNoBounce = [](const EVADE::Fieldspace::TunnellingDir &a,
                           const EVADE::Fieldspace::TunnellingDir &b) {
      return static_cast<EVADE::Fieldspace::Direction>(a) <
             static_cast<EVADE::Fieldspace::Direction>(b);
    };
    const auto glob =
        *std::min_element(result.begin(), result.end(), lessNoBounce);
    storeFSP(glob, "GM");

    return EVADE::CheckStability(mdm);
  }

private:
  using EModel = typename Model::ModelEVADE;
  const std::vector<std::vector<std::string>> fieldSets_;
  EVADE::StationarityConditions<EModel> statConds_;
};

} // namespace ScannerS::Constraints
