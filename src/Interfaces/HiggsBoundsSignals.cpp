#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"

namespace ScannerS::Interfaces::HiggsBoundsSignals {
double tWHratio(double c_HWW, std::complex<double> c_Htt) {
  return (-62. * c_HWW * c_Htt.real() + 45. * pow(c_Htt.real(), 2) +
          34. * pow(c_Htt.imag(), 2) + 33. * pow(c_HWW, 2)) /
         50.;
}
} // namespace ScannerS::Interfaces::HiggsBoundsSignals
