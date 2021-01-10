#pragma once

#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Interfaces/MicrOMEGAs/MicromegasInterface.hpp"
#include <type_traits>

namespace ScannerS::Constraints {

//! functions related to the DarkMatter constraint
namespace DarkMatterDetail {
/**
 * @brief Fit function to the Xenon1T bound on the direct detection cxn.
 *
 * From Fig. 5 of [1805.12562](https://arxiv.org/abs/1805.12562) slightly
 * extrapolated to cover a larger mass range.
 *
 * @param m mass in GeV
 * @return double \f$ 2\sigma \f$ upper bound on the cxn **in pb**
 * @relatedalso ScannerS::Constraints::DarkMatter
 */
double Xenon1TBound(double m);

/**
 * @brief Get the DM-Xe cross section from the DM-nucleon cross sections
 *
 *
 * @param cxn dark matter nucleon cross sections from MicrOMEGAs
 * @return double DM-XE scattering cross section in pb
 * @relatedalso ScannerS::Constraints::DarkMatter
 */
double XenonSICxn(const ScannerS::Interfaces::MicrOMEGAs::DDCxn &cxn);
} // namespace DarkMatterDetail

/**
 * @brief Constraint from dark matter observables.
 *
 * Uses [MicrOMEGAs](https://lapth.cnrs.fr/micromegas/) to calculate DM
 * observables and checks the upper bound on the relic density and bounds from
 * direct detection.
 *
 * @tparam Model a model class with a `Model::micromegasModelName` that
 * corresponds to a folder in `src/Interfaces/MicrOMEGAS/` containing the
 * CalcHEP model files. Additionally a
 * `Model::MOInput(p)->std::map<std::string,double>` function is required, that
 * returns the values of the model parameters using the names in the `vars1.mdl`
 * as keys for a `Model::ParameterPoint` `p`.
 */
template <class Model> class DarkMatter : public Constraint<DarkMatter, Model> {
public:
  static constexpr auto constraintId = "DM"; //!< unique constraint ID

  //! lower bound of the direct detection mass range in GeV
  static constexpr auto ddMassMin = 1;
  //! upper bound of the direct detection mass range in GeV
  static constexpr auto ddMassMax = 2e3;

  //! observed DM density by Planck 2018
  //! [1807.06209](https://arxiv.org/abs/1807.06209)
  static constexpr double omegaC = 0.1200;
  //! uncertainty on the observed DM density by Planck 2018
  //! [1807.06209](https://arxiv.org/abs/1807.06209)
  static constexpr double sdOmegaC = 0.0012;

  /**
   * @brief Constructor that initializes MicrOMEGAs to the correct model
   *
   * **NEVER** create more than one of these objects, MicrOMEGAs can't handle
   * that.
   *
   * @param severity the constraint severity
   * @param relativeDDMassRes approximate relative mass resolution of the direct
   * detection experiment, in multicomponent DM, cross sections are summed if
   * the components are within this mass range
   */
  explicit DarkMatter(Severity severity, double relativeDDMassRes = 0.2)
      : Constraint<DarkMatter, Model>{severity}, relativeDDMassRes_{
                                                     relativeDDMassRes} {
    Interfaces::MicrOMEGAs::SelectModel(Model::micromegasModelName);
  }

  /**
   * @brief Calculates DM observables using MicrOMEGAs and applies the
   * corresponding constraints.
   *
   * Automatically detects whether DM is one or two component. Stores the relic
   * density `omega_c` as well as the spin dependent (`Sd`) and spin independent
   * (`Si`) proton and neutron DM-scattering cross sections in `pb` (names
   * starting with `DD`). Also stores the mass of the DM particle(s) `mDM...`
   * and, for two component DM, the relic density fraction of the second
   * component (`fracDM2`).
   *
   * Requires the relic density not to oversaturate the observed DM density and
   * imposes a bound from direct detection.
   *
   * @param p the parameter point
   * @return is `p` excluded by dark matter observables
   */
  bool Apply(typename Model::ParameterPoint &p) {
    using namespace Interfaces::MicrOMEGAs;
    AssignMOValues(Model::MOInput(p));
    const auto [qDM1, qDM2] = FindDMCandidates();
    if ((qDM1.chargeX3 != 0) || (qDM1.colorDim != 1)) {
      return false;
    } else if ((qDM2.mass > 0) &&
               ((qDM2.chargeX3 != 0) || (qDM2.colorDim != 1))) {
      return false;
    }

    const auto [omega, fracCDM2] = RelicDensity();
    const auto [DD1, DD2] = DDCrossSections();
    p.data.Store("omega_c", omega);
    if (fracCDM2 > 0) {
      p.data.Store("DD1_pSi", DD1.pSi);
      p.data.Store("DD1_pSd", DD1.pSd);
      p.data.Store("DD1_nSi", DD1.nSi);
      p.data.Store("DD1_nSd", DD1.nSd);
      p.data.Store("mDM1", qDM1.mass);
      p.data.Store("DD2_pSi", DD2.pSi);
      p.data.Store("DD2_pSd", DD2.pSd);
      p.data.Store("DD2_nSi", DD2.nSi);
      p.data.Store("DD2_nSd", DD2.nSd);
      p.data.Store("mDM2", qDM2.mass);
      p.data.Store("fracDM2", fracCDM2);
    } else {
      p.data.Store("DD_pSi", DD1.pSi);
      p.data.Store("DD_pSd", DD1.pSd);
      p.data.Store("DD_nSi", DD1.nSi);
      p.data.Store("DD_nSd", DD1.nSd);
      p.data.Store("mDM", qDM1.mass);
    }
    if (omega < 0) // error in micromegas
      return false;

    // simply add the cxns if the particles are close in mass (massRes)
    if (const double avgmass = (qDM1.mass + qDM2.mass) / 2.;
        std::abs(qDM1.mass - qDM2.mass) / avgmass < relativeDDMassRes_) {
      return CheckRelicDensity(omega) &&
             CheckDirectDetection(DD1 * (1 - fracCDM2) + DD2 * fracCDM2, omega,
                                  avgmass);
    } else { // otherwise, check independently
      return CheckRelicDensity(omega) &&
             CheckDirectDetection(DD1, omega * (1 - fracCDM2), qDM1.mass) &&
             CheckDirectDetection(DD2, omega * fracCDM2, qDM2.mass);
    }
  }

private:
  bool CheckRelicDensity(double omega) { return omega < omegaC + 2 * sdOmegaC; }

  bool CheckDirectDetection(const Interfaces::MicrOMEGAs::DDCxn &cxn,
                            double omega, double mass) {
    if (mass < ddMassMin || mass > ddMassMax)
      return true;
    return DarkMatterDetail::XenonSICxn(cxn) * omega / omegaC <
           DarkMatterDetail::Xenon1TBound(mass);
  }

  double relativeDDMassRes_;
};

} // namespace ScannerS::Constraints
