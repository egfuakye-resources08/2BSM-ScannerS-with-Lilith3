#include "ScannerS/Constants.hpp"
#include "ScannerS/Models/C2HDM.hpp"
#include "catch.hpp"
#include "prettyprint.hpp"
#include <cmath>

namespace {
bool BFBold(const std::array<double, 6> &L) {
  if (!(L[0] > 0 && L[1] > 0))
    return false;
  if (!(std::sqrt(L[0] * L[1]) + L[2] > 0))
    return false;
  if (!(std::sqrt(L[0] * L[1]) + L[2] + L[3] >
        std::sqrt(L[4] * L[4] + L[5] * L[5])))
    return false;

  return true;
}

double Unitarityold(const std::array<double, 6> &L) {
  using ScannerS::Constants::pi;
  using std::abs;
  using std::pow;
  using std::sqrt;
  return std::max(
      {abs(3 / 2. * (L[0] + L[1]) +
           sqrt(9 / 4. * pow(L[0] - L[1], 2) + pow(2 * L[2] + L[3], 2))),
       abs(3 / 2. * (L[0] + L[1]) -
           sqrt(9 / 4. * pow(L[0] - L[1], 2) + pow(2 * L[2] + L[3], 2))),
       abs(1 / 2. * (L[0] + L[1]) +
           1 / 2. * sqrt(pow(L[0] - L[1], 2) + 4 * L[3] * L[3])),
       abs(1 / 2. * (L[0] + L[1]) -
           1 / 2. * sqrt(pow(L[0] - L[1], 2) + 4 * L[3] * L[3])),
       abs(1 / 2. * (L[0] + L[1]) +
           1 / 2. *
               sqrt(pow(L[0] - L[1], 2) + 4 * (L[4] * L[4] + L[5] * L[5]))),
       abs(1 / 2. * (L[0] + L[1]) -
           1 / 2. *
               sqrt(pow(L[0] - L[1], 2) + 4 * (L[4] * L[4] + L[5] * L[5]))),
       abs(L[2] + 2 * L[3] + 3 * sqrt(L[4] * L[4] + L[5] * L[5])),
       abs(L[2] + 2 * L[3] - 3 * sqrt(L[4] * L[4] + L[5] * L[5])),
       abs(L[2] + sqrt(L[4] * L[4] + L[5] * L[5])),
       abs(L[2] - sqrt(L[4] * L[4] + L[5] * L[5])), abs(L[2] + L[3]),
       abs(L[2] - L[3])});
}
} // namespace

TEST_CASE("C2HDM BFB and unitarity", "[unit][bfb][c2hdm][unitarity]") {
  std::vector<std::array<double, 6>> testpoints{
      {1.02139, 9.84584, -6.22022, 5.44552, 5.80709, -3.39704},
      {7.55009, 6.67431, 9.54183, -9.72435, 3.62375, 5.02166},
      {6.66685, 8.21315, 7.0416, 2.7642, 4.45829, 8.39388},
      {0.442027, -8.11578, 4.97078, 3.22815, 3.73903, -2.02539},
      {-5.23443, 1.43382, -4.04322, -3.91811, -0.345708, 4.11567},
      {-5.87709, -8.50317, -8.70333, 5.48947, -7.84307, -7.61817},
      {1.69773, 8.94409, -0.604418, 5.88233, 6.07752, -4.48743},
      {7.31615, -8.98566, -1.44969, 6.97432, -2.41353, -5.27853},
      {-7.10629, 6.4093, 8.77301, 1.54124, 5.0162, 9.85419},
      {-8.09307, -3.20029, 4.74302, 0.827871, 7.3875, -1.22088}};

  for (const auto &tp : testpoints) {
    REQUIRE(BFBold(tp) == ScannerS::Models::C2HDM::BFB(tp));
    INFO(tp);
    REQUIRE(Approx(Unitarityold(tp)) ==
            ScannerS::Models::C2HDM::MaxUnitarityEV(tp));
  }
}
