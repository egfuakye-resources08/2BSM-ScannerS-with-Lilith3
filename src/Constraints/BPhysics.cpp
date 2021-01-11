#include "ScannerS/Constraints/BPhysics.hpp"

#include <cmath>

namespace ScannerS::Constraints::BPhysicsDetail {

bool T1Bdmumu(double tbeta, double mHp) {
  return tbeta > 7.641904645416302 + 0.00010838775771032054 * mHp -
                     0.9781661940035135 * std::log(mHp);
}

bool T2Bsmumu(double tbeta, double mHp) {
  return tbeta < 5.305277092184595 + 0.02077277078469879 * mHp;
}

bool T2Bsgam(double tbeta, double mHp) {
  return mHp > 590.4800674010913 - 5.838933724168528 / std::pow(tbeta, 3) +
                   51.863388199691144 / std::pow(tbeta, 2);
}

bool LSBdmumu(double tbeta, double mHp) {
  return tbeta > 7.847284035235539 + 0.00016810529067776107 * mHp -
                     1.018923198977455 * std::log(mHp);
}

bool FBdmumu(double tbeta, double mHp) {
  return tbeta > 7.88241222336676 + 0.0002454822798956865 * mHp -
                     1.0259515796428502 * std::log(mHp);
}

bool FBsgam(double tbeta, double mHp) {
  return mHp > 590.2753776022662 - 5.893554600834598 / std::pow(tbeta, 3) +
                   56.836760951151746 / std::pow(tbeta, 2);
}

} // namespace ScannerS::Constraints::BPhysicsDetail
