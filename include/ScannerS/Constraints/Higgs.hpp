#pragma once

#include "ScannerS/Constants.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp" // IWYU pramga: export
#include <cstddef>
#include <string>

using namespace std::string_literals;

namespace ScannerS::Constraints {

/**
 * @brief Constraint from Higgs searches and measurements.
 *
 * Uses the [HiggsBounds](https://gitlab.com/higgsbounds/higgsbounds) and
 * [HiggsSignals](https://gitlab.com/higgsbounds/higgssignals) libraries.
 *
 * @tparam Model a model class with a `Model::HiggsBoundsInput(p, hbhs)`
 * function that returns a class compatible with
 * Interfaces::HiggsBoundsSignals::HiggsBoundsSignals::RunHBHS() for a
 * `Model::ParameterPoint` `p` and the
 * `Interfaces::HiggsBoundsSignals::HiggsBoundsSignals` object `hbhs`. The model
 * class also needs to provide `Model::nHzero` and `Model::nHplus` members
 * specifying the number of neutral and charged Higgs bosons, as well as
 * `Model::namesHzero` and `Model::namesHplus` arrays of strings containing
 * corresponding names.
 */
template <class Model> class Higgs : public Constraint<Higgs, Model> {
public:
  static constexpr auto constraintId = "Higgs"; //!< unique constraint ID

  //! reference SM Higgs mass
  static constexpr double mhref = 125.09;
  //! reference SM \f$\chi^2\f$ from rate measurements
  static constexpr double chisqMuSM = 84.4372199363;
  //! reference SM \f$\chi^2\f$ from mass measurements
  static constexpr double chisqMassSM = 0.;

  //! Constructor that sets the severity and \f$ \chi^2_\mathrm{crit} \f$. The
  //! chisqCut argument can be set directly in the main function.
  Higgs<Model>(Severity severity, double chisqCut)
      : Constraint<Higgs, Model>{severity}, _chisqCut{chisqCut} {}

  /**
   * @brief Obtains the constraints from Higgs searches and Higgs measurements.
   *
   * Stores `hb_result` and `hb_channel` for the overall HiggsBounds results as
   * well as `hb_Hi_result`, `hb_Hi_channel`, `hb_Hi_obsratio`, and
   * `hb_Hi_ncombined` for each Higgs boson, where `Hi` is replaced by the
   * corresponding name (from `Model::namesHzero` and `Model::namesHplus`).
   *
   * Regarding HiggsSignals results, the mass and rate \f$ \chi^2 \f$
   * contributions are stored as `hs_chisqMass` and `hs_chisqMu`, while
   * \f$\chi^2 - \chi^2_\mathrm{SM} \f$ is stored as `hs_deltaChisq`.
   * Additionally, for each neutral Higgs boson the SM-normalized signal rates
   * in the most important measured channels are stored as `mu_Hi_WW`,
   * `mu_Hi_ZZ`, `mu_Hi_gamgam`, `mu_Hi_tautau`, `mu_Hi_bb`, `mu_Hi_bb_VH`
   * (where Hi is replaced by the entries of `Model::namesHzero`).
   *
   * If needed, additional information can be stored within the
   * `Model::HiggsBoundsInput` function. This is particularly useful for storing
   * the cross sections tabulated within HiggsBounds, as the HiggsBoundsSignals
   * object is available in this function. This can eg be used to store the
   * charged Higgs cross section `x_tHpm`.
   *
   * @param p the parameter point
   * @return is `p` compatible with all Higgs data?
   */
  bool Apply(typename Model::ParameterPoint &p) {
    auto hbin = Model::HiggsBoundsInput(p, hbhs_);
    auto res = hbhs_.RunHBHS(hbin);

    p.data.Store("hb_result", res.result[0]);
    p.data.Store("hb_channel", res.chan[0]);
    for (size_t i = 0; i != HiggsBS::nHzero; ++i) {
      p.data.Store("hb_"s + Model::namesHzero[i] + "_result",
                   res.result[i + 1]);
      p.data.Store("hb_"s + Model::namesHzero[i] + "_channel", res.chan[i + 1]);
      p.data.Store("hb_"s + Model::namesHzero[i] + "_obsratio",
                   res.obsratio[i + 1]);
      p.data.Store("hb_"s + Model::namesHzero[i] + "_ncombined",
                   res.ncombined[i + 1]);
    }
    for (size_t i = 0; i != HiggsBS::nHplus; ++i) {
      p.data.Store("hb_"s + Model::namesHplus[i] + "_result",
                   res.result[i + HiggsBS::nHzero + 1]);
      p.data.Store("hb_"s + Model::namesHplus[i] + "_channel",
                   res.chan[i + HiggsBS::nHzero + 1]);
      p.data.Store("hb_"s + Model::namesHplus[i] + "_obsratio",
                   res.obsratio[i + HiggsBS::nHzero + 1]);
      p.data.Store("hb_"s + Model::namesHplus[i] + "_ncombined",
                   res.ncombined[i + HiggsBS::nHzero + 1]);
    }
    p.data.Store("hs_chisqMu", res.chisqMu);
    p.data.Store("hs_chisqMass", res.chisqMass);
    double deltaChisq = res.chisq - chisqMuSM - chisqMassSM;
    p.data.Store("hs_deltaChisq", deltaChisq);
    for (size_t i = 0; i != HiggsBS::nHzero; ++i) {
      p.data.Store("mu_"s + Model::namesHzero[i] + "_WW", res.mu_WW[i]);
      p.data.Store("mu_"s + Model::namesHzero[i] + "_ZZ", res.mu_ZZ[i]);
      p.data.Store("mu_"s + Model::namesHzero[i] + "_gamgam", res.mu_gaga[i]);
      p.data.Store("mu_"s + Model::namesHzero[i] + "_tautau", res.mu_tautau[i]);
      p.data.Store("mu_"s + Model::namesHzero[i] + "_bb", res.mu_bb[i]);
      p.data.Store("mu_"s + Model::namesHzero[i] + "_bb_VH", res.mu_bb_VH[i]);
    }
    return (res.result[0] == 1) && (deltaChisq < _chisqCut);
  }

private:
  using HiggsBS =
      Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<Model::nHzero,
                                                         Model::nHplus>;
  HiggsBS hbhs_{};
  double _chisqCut;
};

} // namespace ScannerS::Constraints
