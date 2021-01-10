#include "ScannerS/Constants.hpp"
#include "ScannerS/Models/R2HDM.hpp"
#include "catch.hpp"

TEST_CASE("R2HDM couplings", "[unit][r2hdm]") {
  using namespace ScannerS::Models;
  ScannerS::Models::R2HDM::AngleInput in{200,
                                         125,
                                         300,
                                         250,
                                         0.742282,
                                         2.5,
                                         1000,
                                         ScannerS::Models::R2HDM::Yuk::typeI,
                                         ScannerS::Constants::vEW};

  ScannerS::Models::R2HDM::ParameterPoint p1{in};
  in.type = R2HDM::Yuk::typeII;
  ScannerS::Models::R2HDM::ParameterPoint p2{in};
  in.type = R2HDM::Yuk::leptonSpecific;
  ScannerS::Models::R2HDM::ParameterPoint p3{in};
  in.type = R2HDM::Yuk::flipped;
  ScannerS::Models::R2HDM::ParameterPoint p4{in};

  R2HDM::CalcCouplings(p1);
  R2HDM::CalcCouplings(p2);
  R2HDM::CalcCouplings(p3);
  R2HDM::CalcCouplings(p4);

  SECTION("gauge couplings independent of type") {
    REQUIRE(p1.data["c_HlVV"] == p2.data["c_HlVV"]);
    REQUIRE(p2.data["c_HhVV"] == p3.data["c_HhVV"]);

    REQUIRE(p1.data["c_HhAZ"] == p2.data["c_HhAZ"]);
    REQUIRE(p3.data["c_HlAZ"] == p4.data["c_HlAZ"]);

    REQUIRE(p1.data["c_HlHpHm"] == p2.data["c_HlHpHm"]);
    REQUIRE(p3.data["c_HlHpHm"] == p4.data["c_HlHpHm"]);
  }

  SECTION("gauge coupling sum rules") {
    REQUIRE(pow(p1.data["c_HlVV"], 2) + pow(p2.data["c_HhVV"], 2) ==
            Approx(1));
  }
  SECTION("u-type couplings independent of type") {
    REQUIRE(p1.data["c_Hluu_e"] == p2.data["c_Hluu_e"]);
    REQUIRE(p2.data["c_Hhuu_e"] == p3.data["c_Hhuu_e"]);
    REQUIRE(p4.data["c_Auu_o"] == p1.data["c_Auu_o"]);
  }

  SECTION("type dependent fermion couplings") {
    REQUIRE(p1.data["c_Hluu_e"] == p1.data["c_Hlll_e"]);
    REQUIRE(p1.data["c_Hhdd_e"] == p1.data["c_Hhll_e"]);
    REQUIRE(p3.data["c_Auu_o"] == -p3.data["c_Add_o"]);
    REQUIRE(p4.data["c_Hluu_e"] == p4.data["c_Hlll_e"]);

    REQUIRE(p1.data["c_Hluu_e"] == p3.data["c_Hldd_e"]);
    REQUIRE(p1.data["c_Hhuu_e"] == p4.data["c_Hhll_e"]);
    REQUIRE(p2.data["c_Add_o"] == p4.data["c_Add_o"]);
  }

  SECTION("check some values against mathematica") {
    CHECK(p1.alpha == 0.742282);
    CHECK(p2.L[0] == Approx(3.15263));
    CHECK(p3.L[1] == Approx(0.504435));
    CHECK(p4.L[2] == Approx(2.59488));
    CHECK(p1.L[3] == Approx(-0.529492));
    CHECK(p2.L[4] == Approx(-1.43672));

    CHECK(p3.data["c_HlHpHm"] == Approx(-139.375));
    CHECK(p4.data["c_HhHpHm"] == Approx(-719.834));
  }
}
