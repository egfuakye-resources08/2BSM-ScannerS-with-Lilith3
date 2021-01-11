#include "ScannerS/Constraints/Constraint.hpp"

#include <iostream>

namespace ScannerS {
namespace Constraints {
std::istream &operator>>(std::istream &in, Severity &sev) {
  int i;
  in >> i;
  if (i >= -1 && i <= 1) {
    sev = static_cast<ScannerS::Constraints::Severity>(i);
    return in;
  }
  throw std::invalid_argument(
      "Severities must take values in 1 (apply), 0 (ignore), -1 (skip).");
}

std::ostream &operator<<(std::ostream &out, Severity &sev) {
  out << static_cast<int>(sev);
  return out;
}
} // namespace Constraints
} // namespace ScannerS
