#include "ScannerS/Models/CxSMDark.hpp"

#include "ScannerS/Constants.hpp"
#include "catch.hpp"

using ScannerS::Models::CxSMDark;

TEST_CASE("CxSM Dark Generate", "[unit][CxSM]") {
  for (double alpha : {-1., -0.5, 0., 0.5, 1.}) {
    auto in = CxSMDark::AngleInput{125, 200, 300, alpha, 246, 150};
    auto p = CxSMDark::ParameterPoint(in);

    REQUIRE(std::sqrt(
                p.L[0] / 4. * p.v * p.v + p.L[1] / 4. * p.vs * p.vs -
                1 / 4. *
                    std::sqrt(
                        std::pow(p.L[0] * p.v * p.v - p.L[1] * p.vs * p.vs, 2) +
                        std::pow(2 * p.L[2] * p.v * p.vs, 2))) ==
            Approx(p.mHl));
    REQUIRE(std::sqrt(
                p.L[0] / 4. * p.v * p.v + p.L[1] / 4. * p.vs * p.vs +
                1 / 4. *
                    std::sqrt(
                        std::pow(p.L[0] * p.v * p.v - p.L[1] * p.vs * p.vs, 2) +
                        std::pow(2 * p.L[2] * p.v * p.vs, 2))) ==
            Approx(p.mHh));
    REQUIRE(p.alpha == Approx(in.alpha));

    SECTION("Switched input ordering") {
      using ScannerS::Constants::pi;
      auto alpha2 = -std::acos(std::sin(in.alpha));
      auto p2 = CxSMDark::ParameterPoint(
          CxSMDark::AngleInput{in.mHb, in.mHa, in.mHX, alpha2, in.v, in.vs});

      REQUIRE(
          std::sqrt(p2.L[0] / 4. * p2.v * p2.v + p2.L[1] / 4. * p2.vs * p2.vs -
                    1 / 4. *
                        std::sqrt(std::pow(p2.L[0] * p2.v * p2.v -
                                               p2.L[1] * p2.vs * p2.vs,
                                           2) +
                                  std::pow(2 * p2.L[2] * p2.v * p2.vs, 2))) ==
          Approx(p2.mHl));
      REQUIRE(
          std::sqrt(p2.L[0] / 4. * p2.v * p2.v + p2.L[1] / 4. * p2.vs * p2.vs +
                    1 / 4. *
                        std::sqrt(std::pow(p2.L[0] * p2.v * p2.v -
                                               p2.L[1] * p2.vs * p2.vs,
                                           2) +
                                  std::pow(2 * p2.L[2] * p2.v * p2.vs, 2))) ==
          Approx(p2.mHh));
      INFO(alpha2);
      CHECK(p.alpha == Approx(p2.alpha));
    }
  }
}
