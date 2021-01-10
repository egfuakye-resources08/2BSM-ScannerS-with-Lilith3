#include "ScannerS/Constants.hpp"
#include "ScannerS/Models/N2HDMDarkS.hpp"
#include "catch.hpp"

TEST_CASE("N2HDM DarkSinglet Angle Input", "[unit][N2HDMDS]") {
  ScannerS::Models::N2HDMDarkS::AngleInput in;
  in.mHa = 373.4616559;
  in.mHb = 998.3090197;
  in.mA = 572.81670694;
  in.mHp = 235.1848166;
  in.mHD = 488.29581638;
  in.tbeta = 2.0;
  in.alpha = asin(0.99894);
  in.m12sq = 10000.;
  in.L6 = 6;
  in.L7 = 7;
  in.L8 = 8;
  in.type = ScannerS::Models::N2HDMDarkS::Yuk::typeI;
  in.v = ScannerS::Constants::vEW;

  ScannerS::Models::N2HDMDarkS::ParameterPoint p{in};
  REQUIRE(p.mHi[0] == in.mHa);
  REQUIRE(p.mHi[1] == in.mHb);
  REQUIRE(p.mA == in.mA);
  REQUIRE(p.mHp == in.mHp);
  REQUIRE(p.mHD == in.mHD);
  REQUIRE(p.tbeta == in.tbeta);
  REQUIRE(p.alpha == in.alpha);
  REQUIRE(p.R(2, 2) == Approx(1));
  REQUIRE(p.R(0, 2) == Approx(0));
  REQUIRE(p.R(1, 2) == Approx(0));
  REQUIRE(p.R(2, 1) == Approx(0));
  REQUIRE(p.R(2, 0) == Approx(0));
  REQUIRE(p.R(0, 0) == Approx(-sin(p.alpha)));
  REQUIRE(p.R(1, 0) == Approx(cos(p.alpha)));
  REQUIRE(p.R(0, 1) == Approx(cos(p.alpha)));
  REQUIRE(p.R(1, 1) == Approx(sin(p.alpha)));
  REQUIRE(p.R(0, 0) * p.R(0, 1) == Approx(-p.R(1, 0) * p.R(1, 1)));
  CHECK(p.L[0] == Approx(10.).margin(1e-2));
  CHECK(p.L[1] == Approx(20.).margin(0.5));
  CHECK(p.L[2] == Approx(3.).margin(1e-1));
  CHECK(p.L[3] == Approx(4.).margin(1e-3));
  CHECK(p.L[4] == Approx(-5.).margin(1e-3));
  CHECK(p.L[5] == in.L6);
  CHECK(p.L[6] == in.L7);
  CHECK(p.L[7] == in.L8);
  CHECK(p.m12sq == in.m12sq);
  CHECK(p.m11sq == Approx(-89122.8).epsilon(2e-2));
  CHECK(p.m22sq == Approx(-492115).epsilon(3e-2));
  CHECK(p.mssq == Approx(2000).epsilon(1e-2));
}
