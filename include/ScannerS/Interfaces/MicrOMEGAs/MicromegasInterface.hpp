#pragma once

#include <map>
#include <string>
#include <string_view>
#include <utility>

namespace ScannerS {
namespace Interfaces {

//! Namespace of the MicrOMEGAs interface
namespace MicrOMEGAs {

/**
 * Sets the MicrOMEGAs global variables to the model specified by modelName.
 * @param modelName name of the MicrOMEGAs model directory
 */
void SelectModel(std::string_view modelName);

/**
 * Assigns values to the MicrOMEGAs parameters. Any parameters in the
 * `vars1.mdl` can be set. **Make sure that the names match, if any value here
 * does not correspond to to a variable in `vars1.mdl` there will be no error.**
 *
 * @param values map of key, value pairs for the parameters
 */
void AssignMOValues(const std::map<std::string, double> &values);

//! Quantum numbers of a DM candidate
struct QuantumNumbers {
  int spinX2;   //!< DM spin * 2
  int chargeX3; //!< DM electric charge * 3
  int colorDim; //!< dimension of the DM color representation
  double mass;  //!< mass of the DM candidate
};

/**
 * Sorts the dark sector and gets the quantum numbers of the DM candidate(s).
 * Has to be called after AssignMOValues and before the MicrOMEGAs evaluation
 * functions.
 * @return pair of structs with the quantum numbers, if only one dark sector
 * exists the the second element is default constructed
 */
std::pair<QuantumNumbers, QuantumNumbers> FindDMCandidates();

//! relic density
struct Relic {
  double omega_c;  //!< relic density
  double fracCDM2; //!< fraction of DM2 relics
};

/**
 * Calculates the DM relic density.
 * @return value of the relic density
 */
Relic RelicDensity();

//! All direct detection cross sections calculated by MicrOMEGAs in pb.
struct DDCxn {
  double pSi; //!< spin-independent DM p scattering
  double nSi; //!< spin-independent DM n scattering
  double pSd; //!< spin-dependent DM p scattering
  double nSd; //!< spin-dependent DM n scattering
};

//! Elementwise sum of two DDCxn
DDCxn operator+(DDCxn a, const DDCxn &b);
//! Elementwise multiplication of a DDCxn with a number
DDCxn operator*(DDCxn cxn, double num);
//! @copydoc operator*(DDCxn,double)
inline DDCxn operator*(double num, DDCxn cxn) { return cxn * num; }

/**
 * Calculates the direct DM detection cross sections.
 * @return DD cross sections
 */
std::pair<DDCxn, DDCxn> DDCrossSections();

} // namespace MicrOMEGAs
} // namespace Interfaces
} // namespace ScannerS
