#include "ScannerS/Constants.hpp"
#include "ScannerS/Models/N2HDMDarkD.hpp"
#include "catch.hpp"

TEST_CASE("N2HDM DarkDoublet Angle Input", "[unit][N2HDMDD]") {
  ScannerS::Models::N2HDMDarkD::AngleInput in;
  in.mHa = 287.01;
  in.mHb = 1422.63;
  in.mAD = 1172.07;
  in.mHDp = 1049.25;
  in.mHD = 1034.71;
  in.alpha = asin(-0.519447);
  in.m22sq = 10000.;
  in.L2 = 20;
  in.L8 = 8;
  in.vs = 500;
  in.v = ScannerS::Constants::vEW;

  ScannerS::Models::N2HDMDarkD::ParameterPoint p{in};
  REQUIRE(p.mHi[0] == in.mHa);
  REQUIRE(p.mHi[1] == in.mHb);
  REQUIRE(p.mAD == in.mAD);
  REQUIRE(p.mHDp == in.mHDp);
  REQUIRE(p.mHD == in.mHD);
  REQUIRE(p.alpha == in.alpha);
  REQUIRE(p.R(2, 1) == Approx(1));
  REQUIRE(p.R(2, 0) == Approx(0));
  REQUIRE(p.R(2, 2) == Approx(0));
  REQUIRE(p.R(0, 1) == Approx(0));
  REQUIRE(p.R(1, 1) == Approx(0));
  REQUIRE(p.R(0, 0) == Approx(cos(p.alpha)));
  REQUIRE(p.R(0, 2) == Approx(sin(p.alpha)));
  REQUIRE(p.R(1, 0) == Approx(-sin(p.alpha)));
  REQUIRE(p.R(1, 2) == Approx(cos(p.alpha)));
  CHECK(p.L[0] == Approx(10.).margin(1e-3));
  CHECK(p.L[1] == in.L2);
  CHECK(p.L[2] == Approx(3.).margin(1e-3));
  CHECK(p.L[3] == Approx(4.).margin(1e-3));
  CHECK(p.L[4] == Approx(-5.).margin(1e-3));
  CHECK(p.L[5] == Approx(6));
  CHECK(p.L[6] == Approx(7));
  CHECK(p.L[7] == in.L8);
  CHECK(p.m11sq == Approx(-1.17812e6));
  CHECK(p.m22sq == in.m22sq);
  CHECK(p.mssq == Approx(-962183).epsilon(1e-3));
}
