#include "ScannerS/Constraints/DarkMatter.hpp"
#include "ScannerS/Constants.hpp"
#include <cmath>

namespace ScannerS::Constraints::DarkMatterDetail {

double Xenon1TBound(double m) {
  using std::pow;
  double logM = std::log10(m);
  double exponent = -47.671797276126775 + 43.43457037808014 / pow(logM, 5) -
                    200.86181549108545 / pow(logM, 4) +
                    364.0472284553867 / pow(logM, 3) -
                    323.67181156510736 / pow(logM, 2) +
                    151.11688814933748 / logM + 4.337753083118923 * logM;
  return pow(10, exponent);
}

double XenonSICxn(const ScannerS::Interfaces::MicrOMEGAs::DDCxn &cxn) {
  using ScannerS::Constants::aNRXe;
  using ScannerS::Constants::aWeightXe;
  return (cxn.pSi * aNRXe + cxn.nSi * (aWeightXe - aNRXe)) / aWeightXe;
}
} // namespace ScannerS::Constraints::DarkMatterDetail
