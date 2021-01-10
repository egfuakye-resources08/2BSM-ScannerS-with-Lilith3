#include "ScannerS/Utilities.hpp"
#include "catch.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace {
const std::vector<std::array<double, 3>> inputangles{
    {-1.39744, -0.0399273, 0.418412},
    {3.07205, 2.53568, -0.997093},
    {2.7355, -1.3225, -2.17389},
    {-1.07772, -2.53117, 0.480547},
    {2.66136, -1.73409, 2.67975},
    {1.19273, -0.281784, 0.836052},
    {2.05455, -2.07747, 2.22036},
    {2.14972, 2.50133, -0.784693},
    {-1.43949, 0.467904, 1.45663},
    {-1.84194, -0.97857, -2.08504},
    {0.0298794, -3.0619, -0.768827},
    {3.05189, -1.8305, -1.10864},
    {-1.31553, 0.0114572, 0.791876},
    {2.55236, -2.79192, -2.02981},
    {-1.91567, 1.67234, -3.11597},
    {1.07362, -1.16747, 0.878817},
    {0.660272, -0.254116, 2.27672},
    {1.17447, -0.403116, -2.20701},
    {0.692055, -1.86222, -1.01534},
    {1.41734, 2.24005, -2.99207},
    {-1.49311, 2.33848, -2.03513},
    {1.23527, -1.33722, -1.78911},
    {1.16894, -1.37506, -0.511868},
    {-2.6966, -1.86174, 1.50837},
    {2.7983, 1.96205, -3.0557},
    {-2.18467, 2.49425, 1.82565},
    {3.07288, -2.17112, 1.12696},
    {1.72175, -2.33517, 2.18906},
    {-2.02867, 0.746004, 0.828995},
    {1.87899, -1.43797, 2.81886},
    {2.84921, 1.54627, 2.27063},
    {-0.187255, -0.379304, -0.266117},
    {1.14245, -0.744216, 2.58728},
    {1.55127, -2.35823, -1.93536},
    {3.09275, -0.91538, -0.889537},
    {2.80077, -2.12026, -0.292163},
    {2.09736, -2.14369, 2.72208},
    {0.430302, -0.300579, -1.83354},
    {0.750968, -0.0855564, -1.57858},
    {-0.53037, 2.87859, -2.34959},
    {0.261274, 1.09519, -0.445055},
    {0.99706, -0.0339228, -1.68129},
    {0.45714, 0.113557, 1.7963},
    {2.61055, -2.59394, -0.117292},
    {-2.62836, 2.60311, 1.74076},
    {-0.774518, -2.42461, -1.93231},
    {0.0622336, 2.04706, 1.53827},
    {0.381179, -2.27175, -0.0218464},
    {0.133141, 0.0813273, 0.777941},
    {0.631506, 2.72164, -0.580822},
    {1.55, 0.8, 1.55}};

const std::vector<std::array<double, 3>> permutations{
    {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};

} // namespace

TEST_CASE("Index Sort", "[unit][utils]") {
  for (const auto &x : inputangles) {
    auto indices = ScannerS::Utilities::IndexSort(x);
    std::array<double, 3> test{x};
    std::sort(test.begin(), test.end());
    for (size_t i = 0; i != x.size(); ++i) {
      REQUIRE(x[indices[i]] == Approx(test[i]));
    }
  }
}

TEST_CASE("Mixmat Parametrization", "[unit][utils]") {
  const double a1 = 0.3;
  const double a2 = -1.1;
  const double a3 = 0.85;
  Eigen::Matrix3d ref;
  ref << 0.433337, 0.134047, -0.891207, 0.444604, 0.828371, 0.340778, 0.78393,
      -0.543906, 0.299366;
  ScannerS::Utilities::MixMat3d(a1, a2, a3).isApprox(ref);
}

TEST_CASE("Ordering Mixing Matrix", "[unit][utils]") {
  for (const auto &ia : inputangles) {
    Eigen::Matrix3d test{Eigen::AngleAxisd(-ia[2], Eigen::Vector3d::UnitX()) *
                         Eigen::AngleAxisd(ia[1], Eigen::Vector3d::UnitY()) *
                         Eigen::AngleAxisd(-ia[0], Eigen::Vector3d::UnitZ())};

    for (const auto &perm : permutations) {
      test =
          Eigen::PermutationMatrix<3>{
              Eigen::Map<Eigen::Vector3i>{
                  ScannerS::Utilities::IndexSort(perm).data()}}
              .transpose() *
          test;
      ScannerS::Utilities::MixMatNormalForm3d(test);
      auto angles = ScannerS::Utilities::MixMatAngles3d(test);
      INFO("ia1 = " + std::to_string(ia[0]));
      INFO("ia2 = " + std::to_string(ia[2]));
      INFO("ia3 = " + std::to_string(ia[3]));
      INFO("perm = (" + std::to_string(perm[0]) + std::to_string(perm[1]) +
           std::to_string(perm[2]) + ")");
      REQUIRE(test.isApprox(Eigen::Matrix3d{
          Eigen::AngleAxisd(-angles[2], Eigen::Vector3d::UnitX()) *
          Eigen::AngleAxisd(angles[1], Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(-angles[0], Eigen::Vector3d::UnitZ())}));
    }
  }
}
