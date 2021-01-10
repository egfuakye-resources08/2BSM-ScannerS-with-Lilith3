#include "ScannerS/Models/R2HDM.hpp"

#include "ScannerS/Constants.hpp"
#include "catch.hpp"

TEST_CASE("R2HDM Generation", "[R2HDM][unit]") {
  using namespace ScannerS::Models;
  SECTION("Correct ordering") {
    R2HDM::AngleInput in{125.09,
                         100,
                         200,
                         200,
                         1.24031081,
                         2.5,
                         2000,
                         R2HDM::Yuk::typeI,
                         ScannerS::Constants::vEW};
    R2HDM::ParameterPoint p{in};
    REQUIRE(p.mHl == 100);
    REQUIRE(p.mHh == 125.09);
    REQUIRE(p.mA == 200);
    REQUIRE(p.mHp == 200);
    REQUIRE(p.alpha == 1.24031081);
    REQUIRE(p.tbeta == 2.5);
    REQUIRE(p.m12sq == 2000);

    CHECK(p.L[0] == Approx(6.69060770e-01));
    CHECK(p.L[1] == Approx(2.72715608e-01));
    CHECK(p.L[2] == Approx(1.30684681e+00));
    CHECK(p.L[3] == Approx(-5.64127725e-01));
    CHECK(p.L[4] == Approx(-5.64127725e-01));

    CHECK(p.m11sq == Approx(-2464.1668706409));
    CHECK(p.m22sq == Approx(-7073.0991108036));
  }
  SECTION("Inverted Ordering") {
    R2HDM::AngleInput in{100,
                         125.09,
                         200,
                         200,
                         -0.330486,
                         2.5,
                         2000,
                         R2HDM::Yuk::typeI,
                         ScannerS::Constants::vEW};
    R2HDM::ParameterPoint p{in};
    REQUIRE(p.mHl == 100);
    REQUIRE(p.mHh == 125.09);
    REQUIRE(p.mA == 200);
    REQUIRE(p.mHp == 200);
    REQUIRE(p.alpha == Approx(1.24031081));
    REQUIRE(p.tbeta == 2.5);
    REQUIRE(p.m12sq == 2000);

    CHECK(p.L[0] == Approx(6.69060770e-01));
    CHECK(p.L[1] == Approx(2.72715608e-01));
    CHECK(p.L[2] == Approx(1.30684681e+00));
    CHECK(p.L[3] == Approx(-5.64127725e-01));
    CHECK(p.L[4] == Approx(-5.64127725e-01));

    CHECK(p.m11sq == Approx(-2464.1668706409));
    CHECK(p.m22sq == Approx(-7073.0991108036));
  }

  SECTION("Physical Input") {
    R2HDM::PhysicalInput in{125.09,
                            100,
                            200,
                            200,
                            -0.05,
                            2.5,
                            2000,
                            R2HDM::Yuk::typeI,
                            ScannerS::Constants::vEW};
    R2HDM::ParameterPoint p{in};
    REQUIRE(p.mHl == 100);
    REQUIRE(p.mHh == 125.09);
    REQUIRE(p.mA == 200);
    REQUIRE(p.mHp == 200);
    REQUIRE(p.alpha == Approx(1.24031081));
    REQUIRE(p.tbeta == 2.5);
    REQUIRE(p.m12sq == 2000);

    CHECK(p.L[0] == Approx(6.69060770e-01));
    CHECK(p.L[1] == Approx(2.72715608e-01));
    CHECK(p.L[2] == Approx(1.30684681e+00));
    CHECK(p.L[3] == Approx(-5.64127725e-01));
    CHECK(p.L[4] == Approx(-5.64127725e-01));

    CHECK(p.m11sq == Approx(-2464.1668706409));
    CHECK(p.m22sq == Approx(-7073.0991108036));

    CHECK(std::sin(std::atan(p.tbeta) - p.alpha) == Approx(-0.05));
  }

  SECTION("Physical Inverted Ordering") {
    R2HDM::PhysicalInput in{100,
                            125.09,
                            200,
                            200,
                            9.98749218e-01,
                            2.5,
                            2000,
                            R2HDM::Yuk::typeI,
                            ScannerS::Constants::vEW};
    R2HDM::ParameterPoint p{in};
    REQUIRE(p.mHl == 100);
    REQUIRE(p.mHh == 125.09);
    REQUIRE(p.mA == 200);
    REQUIRE(p.mHp == 200);
    REQUIRE(p.alpha == Approx(1.24031081));
    REQUIRE(p.tbeta == 2.5);
    REQUIRE(p.m12sq == 2000);

    CHECK(p.L[0] == Approx(6.69060770e-01));
    CHECK(p.L[1] == Approx(2.72715608e-01));
    CHECK(p.L[2] == Approx(1.30684681e+00));
    CHECK(p.L[3] == Approx(-5.64127725e-01));
    CHECK(p.L[4] == Approx(-5.64127725e-01));

    CHECK(p.m11sq == Approx(-2464.1668706409));
    CHECK(p.m22sq == Approx(-7073.0991108036));
    CHECK(std::sin(std::atan(p.tbeta) - p.alpha) == Approx(-0.05));
  }
}
